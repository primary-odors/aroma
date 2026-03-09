<?php

chdir(__DIR__);
include_once("protutils.php");

function lencmp($a, $b)
{
    $ca = is_array($a) ? count($a) : strlen($a);
    $cb = is_array($b) ? count($b) : strlen($b);

    if ($ca == $cb) return 0;
    return ($ca < $cb) ? 1 : -1;
}

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
    // shuffle($keys);

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
    if (isset($sequences["TAAR1"])) $outgroup = "TAAR1";
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

function find_bt_of_root($clade)
{
    $cursor = false;
    foreach($clade as $bt)
    {
        if (!$cursor) $cursor = $bt;
        else
        {
            for ($i=strlen($cursor); $i; $i--)
            {
                if (substr($cursor, 0, $i) == substr($bt, 0, $i))
                {
                    $cursor = substr($cursor, 0, $i);
                    break;
                }
            }
        }
    }
    return $cursor;
}

$ali = [];
$c = file_get_contents("sequences_aligned.txt");
$ells = [];
foreach (explode("\n", $c) as $ln)
{
    $id = trim(substr($ln, 0, 7));
    $sq = substr($ln, 8);

    if (isset($prots[$id]))
    {
        $ali[$id] = $sq;
        // foreach ($ells as $l) $ali[$id] .= substr($sq, $l, 1);
    }
    else if (preg_match("/^[LM\\s]+$/", $ln))
    {
        $ells = [];
        $j=0;
        while ($j = strpos($sq, 'L', $j))
        {
            $ells[] = $j;
            $j++;
        }
    }
}

if (!count($ells)) die("No ells.\n");
// print_r($ells); exit;

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
uasort($consensus, "lencmp");

foreach ($subseqs as $fam => $subs)
{
    foreach ($subs as $sub => $seqs)
    {
        $subsensus[$fam][$sub] = generate_consensus($seqs);
    }
    uasort($subsensus[$fam], "lencmp");
}

$tnodes = [];
$cladeseqs =
[
    "ClassI" => generate_consensus([ $consensus["OR51"], $consensus["OR52"], $consensus["OR56"] ]),
    "OR1/3/7" => generate_consensus([ $consensus["OR1"], $consensus["OR3"], $consensus["OR7"] ]),
    "OR2/13" => generate_consensus([ $consensus["OR2"], $consensus["OR13"] ]),
    "OR4/12" => generate_consensus([ $consensus["OR4"], $consensus["OR12"] ]),
    "OR5/8/9" => generate_consensus([ $consensus["OR5"], $consensus["OR8"], $consensus["OR9"] ]),
    "OR6/10/11" => generate_consensus([ $consensus["OR6"], $consensus["OR10"], $consensus["OR11"] ]),
    "OR14" => $consensus["OR14"],
];

$cladetree = generate_btree($cladeseqs, "00");
uasort($cladetree, "lencmp");
print_r($cladetree);

foreach ($cladetree as $k => $bt) $tnodes["*$bt"] = $k;

$famtree = generate_btree(["OR51" => $consensus["OR51"], "OR52" => $consensus["OR52"], "OR56" => $consensus["OR56"]], $cladetree["ClassI"]);
$famtree = array_merge($famtree, generate_btree(["OR1" => $consensus["OR1"], "OR3" => $consensus["OR3"], "OR7" => $consensus["OR7"]], $cladetree["OR1/3/7"]));
$famtree = array_merge($famtree, generate_btree(["OR2" => $consensus["OR2"], "OR13" => $consensus["OR13"]], $cladetree["OR2/13"]));
$famtree = array_merge($famtree, generate_btree(["OR4" => $consensus["OR4"], "OR12" => $consensus["OR12"]], $cladetree["OR4/12"]));
$famtree = array_merge($famtree, generate_btree(["OR5" => $consensus["OR5"], "OR8" => $consensus["OR8"], "OR9" => $consensus["OR9"]], $cladetree["OR5/8/9"]));
$famtree = array_merge($famtree, generate_btree(["OR6" => $consensus["OR6"], "OR10" => $consensus["OR10"], "OR11" => $consensus["OR11"]], $cladetree["OR6/10/11"]));
$famtree = array_merge($famtree, generate_btree(["OR14" => $consensus["OR14"]], $cladetree["OR14"]));
$famtree = array_merge($famtree, generate_btree(["TAAR" => $consensus["TAAR"]], "01"));
$famtree = array_merge($famtree, generate_btree(["VN1R" => $consensus["VN1R"]], "1"));

uasort($famtree, "lencmp");
print_r($famtree);

foreach ($famtree as $k => $bt) $tnodes["$bt"] = $k;

$subtree = [];

foreach ($subseqs as $fam => $subs)
{
    $subtree[$fam] = generate_btree($subsensus[$fam], $famtree[$fam]);
    uasort($subtree[$fam], "lencmp");
    // print_r($subtree);    exit;
}
uasort($subtree, "lencmp");
print_r($subtree);

foreach ($subtree as $fam => $subs) foreach ($subs as $sub => $bt) $tnodes["$bt"] = $fam.$sub;
$tnodes['*'.find_bt_of_root([$subtree["OR2"]["M"], $subtree["OR2"]["T"], $subtree["OR2"]["V"]])] = "OR2M/T/V";
$tnodes['*'.find_bt_of_root([$subtree["OR5"]["A"], $subtree["OR5"]["AN"]])] = "OR5A/AN";
$tnodes['*'.find_bt_of_root([$subtree["OR6"]["A"], $subtree["OR6"]["B"], $subtree["OR6"]["P"], $subtree["OR6"]["Y"]])] = "Schiff6";
$tnodes['*'.find_bt_of_root([$subtree["OR10"]["D"], $subtree["OR10"]["G"], $subtree["OR10"]["S"]])] = "OR10D/G/S";

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
        $ltree = generate_btree($sequences, $fb);
        $rcptree = array_merge($rcptree, $ltree);
    }
}
// print_r($rcptree);

foreach ($rcptree as $protid => $bt) if (isset($prots[$protid])) $prots[$protid]["btree"] = $bt;
natsort($tnodes);

file_put_contents("../data/receptor.json", json_encode_pretty($prots));
file_put_contents("../data/tree_nodes.json", json_encode_pretty($tnodes));

