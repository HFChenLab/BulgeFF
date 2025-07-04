def atom_to_pdb(atom_info):
    pdb_lines = []
    for atom in atom_info:
        pdb_lines.append(
            f"ATOM  {atom['atom_number']:5d} {atom['atom_name']:<4} {atom['residue_name']:<3} {atom['chain_id']}{atom['residue_number']:4d}    {atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}{atom['occupancy']:6.2f}{atom['bfactor']:6.2f}          {atom['element']:<2}"
        )
    return "\n".join(pdb_lines)