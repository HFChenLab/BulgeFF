def write_plumed(plumed_file, eta_list, theta_list, eta_func_list, theta_func_list, bulge_res_list, bulge_resi_list):
    with open(plumed_file, "w") as f:
        f.write(f"MOLINFO STRUCTURE=reference.pdb MOLTYPE=rna\n")
        for eta, theta, eta_func, theta_func, bulge_res, bulge_resi in zip(eta_list, theta_list, eta_func_list, theta_func_list, bulge_res_list, bulge_resi_list):
            identifier = f"{bulge_res}{bulge_resi}"
            f.write(f"\n# Bulge {identifier}\n")
            f.write(f"eta_{identifier}: TORSION ATOMS={eta}\n")
            f.write(f"theta_{identifier}: TORSION ATOMS={theta}\n")
            f.write(f"bias_eta_{identifier}: CUSTOM ARG=eta_{identifier} FUNC={eta_func} PERIODIC=NO\n")
            f.write(f"bias_theta_{identifier}: CUSTOM ARG=theta_{identifier} FUNC={theta_func} PERIODIC=NO\n")
            f.write(f"bias_e_{identifier}: BIASVALUE ARG=bias_eta_{identifier}\n")
            f.write(f"bias_t_{identifier}: BIASVALUE ARG=bias_theta_{identifier}\n")
            f.write(f"PRINT ARG=eta_{identifier},bias_e_{identifier}.bias FILE=eta_{identifier}.dat STRIDE=5000\n")
            f.write(f"PRINT ARG=theta_{identifier},bias_t_{identifier}.bias FILE=theta_{identifier}.dat STRIDE=5000\n")
    f.close()