<?php

chdir(__DIR__);
include_once("protutils.php");

function generate_consensus($aligneds)
{
    $bins = [];
    foreach ($aligneds as $seq)
    {
        for ($i=0; $i<strlen($seq); $i++)
        {
            if (!isset($bins[$i])) $bins[$i] = [];
            $c = substr($seq, $i, 1);
            if ($c < 'A' || $c > 'V') continue;
            if (!isset($bins[$i][$c])) $bins[$i][$c] = 1;
            else $bins[$i][$c]++;
        }
    }

    $consensus = "";
    foreach ($bins as $i => $b)
    {
        if (!count($b)) $consensus .= "-";
        else
        {
            arsort($b);
            $consensus .= array_keys($b)[0];
        }
    }

    return $consensus;
}

function generate_btree($sequences, $prefix = "")
{
    $cs = count($sequences);
    if ($cs == 1) return [array_keys($sequences)[0] => $prefix];
    else if ($cs == 2)
        return [array_keys($sequences)[0] => $prefix.'0', array_keys($sequences)[1] => $prefix.'1'];

    $keys = array_keys($sequences);
    shuffle($keys);

    $consensus = generate_consensus($sequences);
    $outgroup = false;
    $outdist = 1e9;
    foreach ($sequences as $k => $seq)
    {
        $dist = similar_text($seq, $consensus);
        if ($dist < $outdist)
        {
            $outdist = $dist;
            $outgroup = $k;
        }
    }
    if (isset($sequences["VN1R"])) $outgroup = "VN1R";
    // echo "OUTGROUP: $outgroup\n";

    $fasta = "";
    foreach($keys as $fam) $fasta .= make_fasta($fam, $sequences[$fam]);

    file_put_contents("../tmp/tmp.fasta", $fasta);
    $og = count($sequences)-1;
    $cmd = "clustalw -infile=../tmp/tmp.fasta -tree -seed=8647 -outputtree=nexus";
    // echo "$cmd\n";
    exec($cmd);
    $roottree = <<<opfisehciet
from Bio import Phylo
from io import StringIO

tree = Phylo.read("../tmp/tmp.tre", "nexus")
outgroup_name = "$outgroup"
outgroup_clade = tree.find_clades(outgroup_name)

try:
    outgroup = next(outgroup_clade)
    tree.root_with_outgroup(outgroup)
    print(f"Tree successfully rerooted with {outgroup_name} as the outgroup.")

except StopIteration:
    print(f"Error: Outgroup '{outgroup_name}' not found in the tree.")

Phylo.write(tree, "../tmp/rooted.tre", "nexus")
opfisehciet;
    file_put_contents("roottree.py", $roottree);
    $cmd = "python3 roottree.py";
    // echo "$cmd\n";
    exec($cmd);

    $btree = [];
    $c = file_get_contents("../tmp/rooted.tre");
    foreach (explode("\n", $c) as $ln)
    {
        if (substr(trim($ln), 0, 5) == "Tree ")
        {
            $treedat = explode('=', $ln, 2)[1];
            $nodists = preg_replace("/:[0-9.-]+/", "", $treedat);
            // echo "$nodists\n";

            $cursor = $prefix;
            while ($nodists)
            {
                $c = substr($nodists, 0, 1);
                if ($c == '(')
                {
                    $cursor .= '0';
                    $nodists = substr($nodists, 1);
                    // echo "$cursor\t\t$nodists\n";
                }
                else if ($c == ',')
                {
                    $motherfucker = 1+intval(substr($cursor, -1));
                    if ($motherfucker > 1)
                    {
                        $l = strlen($cursor)-1;
                        foreach ($btree as $k => $v) $btree[$k] = substr($v, 0, $l).'0'.substr($v, $l);
                        $motherfucker--;
                    }
                    $cursor = substr($cursor, 0, -1).$motherfucker;
                    $nodists = substr($nodists, 1);
                    // echo "$cursor\t\t$nodists\n";
                }
                else if ($c == ')')
                {
                    $cursor = substr($cursor, 0, -1);
                    $nodists = substr($nodists, 1);
                    // echo "$cursor\t\t$nodists\n";
                }
                else if ($c == ';')
                {
                    break;
                }
                else
                {
                    $m = [];
                    preg_match("/^[^();,]+/", $nodists, $m);
                    $name = $m[0];
                    $btree[$name] = $cursor;
                    $nodists = substr($nodists, strlen($name));
                }
            }

            break;
        }
    }

    // print_r($btree);
    return $btree;
}

$ali = [];
$c = file_get_contents("sequences_aligned.txt");
foreach (explode("\n", $c) as $ln)
{
    $id = trim(substr($ln, 0, 7));
    $sq = substr($ln, 8);

    if (isset($prots[$id])) $ali[$id] = $sq;
}

$consensus = [];
$subsensus = [];
$famseqs = [];
$subseqs = [];
foreach ($prots as $protid => $p)
{
    if (!isset($ali[$protid])) continue;
    $fam = family_from_protid($protid);
    if (!isset($famseqs[$fam])) $famseqs[$fam] = [];
    $famseqs[$fam][$protid] = $ali[$protid];

    $sub = subfamily_from_protid($protid);
    if (!isset($subseqs[$fam][$sub])) $subseqs[$fam][$sub] = [];
    $subseqs[$fam][$sub][$protid] = $ali[$protid];
}

foreach ($famseqs as $fam => $seqs)
{
    $consensus[$fam] = generate_consensus($seqs);
}

foreach ($subseqs as $fam => $subs)
{
    foreach ($subs as $sub => $seqs)
    {
        $subsensus[$fam][$sub] = generate_consensus($seqs);
    }
}

$famtree = generate_btree($consensus);
print_r($famtree);

$subtree = [];

foreach ($subseqs as $fam => $subs)
{
    $subtree[$fam] = generate_btree($subsensus[$fam], $famtree[$fam]);
    // print_r($subtree);    exit;
}
print_r($subtree);

$rcptree = [];
foreach ($subtree as $fam => $subs)
{
    foreach ($subs as $sub => $fb)
    {
        $sequences = [];
        foreach (array_keys($prots) as $protid)
        {
            if ($fam != family_from_protid($protid)) continue;
            if ($sub != subfamily_from_protid($protid)) continue;
            $sequences[$protid] = $ali[$protid];
        }
        // $sequences["VN1R"] = "DQNTMINDMEGIUSTAGARPAGESEQVENCEDESIGNEDTQALWAYSCARRYTHERQLEQFRQQTNQDE";
        $ltree = generate_btree($sequences, $fb);
        $rcptree = array_merge($rcptree, $ltree);
    }
}
// print_r($rcptree);

foreach ($rcptree as $protid => $bt) if (isset($prots[$protid])) $prots[$protid]["btree"] = $bt;

$fp = fopen("../data/receptor.json", "w");
if ($fp)
{
	fwrite($fp, json_encode_pretty($prots));
	fclose($fp);
}

