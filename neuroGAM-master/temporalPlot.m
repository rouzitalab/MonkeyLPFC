l_sac_avg = [0.008077845, 0.034145954, 0.030661456, 0.048194358];
l_rul_avg = [0.010049469, 0.035681501, 0.032692592, 0.045730481];
l_sac_std = [0.010321519	0.029063508	0.03321726	0.046120333];
l_rul_std = [0.011305933	0.027705481	0.032249659	0.041007726];
l_sac_var = [0.000106534	0.000844688	0.001103386	0.002127085];
l_rul_var = [0.000127824	0.000767594	0.001040041	0.001681634];
ul_sac_avg = [0.013387796	0.040270984	0.031806308	0.051736378];
ul_rul_avg = [0.014545463	0.040844903	0.030491861	0.037702981];
ul_sac_std = [0.024587218	0.035104624	0.032241261	0.042738367];
ul_rul_std = [0.023831716	0.03345331	0.027722307	0.033734355];
ul_sac_var = [0.000604531	0.001232335	0.001039499	0.001826568];
ul_rul_var = [0.000567951	0.001119124	0.000768526	0.001138007];
errorbar(l_sac_avg, l_sac_var)
errorbar(l_rul_avg, l_rul_var)