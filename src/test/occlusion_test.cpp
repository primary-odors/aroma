
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include "../classes/protein.h"

using namespace std;

int main(int argc, char** argv)
{
    Protein p("RECEPTOR");
    Molecule m("LIGAND");
    FILE *fp, *lfp;
    char buffer[1024];
    int i, j;
    char c;

    if (argc < 2)
    {
        cerr << "A dock file is required." << endl << flush;
        return -1;
    }
    std::string dfname = argv[1];

    fp = fopen(dfname.c_str(), "rb");
    if (!fp)
    {
        cerr << "Failed to open " << dfname << " for reading." << endl << flush;
        return -1;
    }

    while (fgets(buffer, 1022, fp))
    {
        char* newline = strchr(buffer, '\n');
        if (newline) *newline = 0;
        newline = strchr(buffer, '\r');
        if (newline) *newline = 0;

        c = buffer[8];
        if (c == ':')
        {
            buffer[8] = 0;
            if (!strcmp(buffer, "PDB file"))
            {
                // cout << "Protein is " << &buffer[10] << endl;
                lfp = fopen(&buffer[10], "rb");
                if (lfp)
                {
                    p.load_pdb(lfp);
                    fclose(lfp);
                    continue;
                }
                else cerr << "Failed to open " << &buffer[10] << " for reading." << endl;
            }
            buffer[8] = c;
        }

        c = buffer[6];
        buffer[6] = 0;
        if (c == ':')
        {
            if (!strcmp(buffer, "Ligand"))
            {
                lfp = fopen(&buffer[8], "rb");
                if (lfp)
                {
                    char sdfdat[16384];
                    fread(sdfdat, 1, 16382, lfp);
                    fclose(lfp);
                    m.from_sdf(sdfdat);
                    continue;
                }
                else cerr << "Failed to open " << &buffer[8] << " for reading." << endl;
            }
        }
        else if (!strcmp(buffer, "HETATM"))
        {
            char aname[8];
            j=0;
            for (i=12; i<16; i++) if (buffer[i] > ' ') aname[j++] = buffer[i];
            aname[j] = 0;
            Atom *a = m.get_atom(aname);
            if (a)
            {
                buffer[54] = 0;
                float z = atof(&buffer[46]);
                buffer[46] = 0;
                float y = atof(&buffer[38]);
                buffer[38] = 0;
                float x = atof(&buffer[30]);
                a->move(Point(x,y,z));
            }
            continue;
        }
        else if (!strcmp(buffer, "ATOM  "))
        {
            char aname[8];
            j=0;
            for (i=12; i<16; i++) if (buffer[i] > ' ') aname[j++] = buffer[i];
            aname[j] = 0;
            buffer[26] = 0;
            int resno = atoi(&buffer[22]);
            if (resno)
            {
                AminoAcid *aa = p.get_residue(resno);
                if (aa)
                {
                    Atom *a = aa->get_atom(aname);
                    if (a)
                    {
                        buffer[54] = 0;
                        float z = atof(&buffer[46]);
                        buffer[46] = 0;
                        float y = atof(&buffer[38]);
                        buffer[38] = 0;
                        float x = atof(&buffer[30]);
                        a->move(Point(x,y,z));
                    }
                }
            }
            continue;
        }
        buffer[6] = c;

        if (buffer[0] == 'E' && buffer[1] == 'N' && buffer[2] == 'D') break;
        if (feof(fp)) break;
    }
    fclose(fp);

    if (!p.get_seq_length())
    {
        cout << "No residues loaded." << endl << flush;
        return 0;
    }
    if (!m.get_atom_count())
    {
        cout << "No ligand atoms loaded." << endl << flush;
        return 0;
    }

    Molecule* mols[SPHREACH_MAX];
    int sphres = p.get_residues_can_clash_ligand((AminoAcid**)mols, &m, m.get_barycenter(), Point(10, 10, 10));
    float o = m.occlusion(mols);

    cout << "Occlusion: " << o << endl;
    return 0;
}