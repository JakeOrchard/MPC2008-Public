def mpctable(data, table_dict, mpcset, simultitle, savename):

    GEset = data[list(mpcset)[0]]['MPCs'].keys()

    f = open('../output/' + savename + '.tex', 'w')
    f.write('\\begin{table}[htbp]\n')
    f.write('\\caption{General Equilibrium Marginal Propensity to Consume: ' + simultitle + '}\n') 
    f.write('\\begin{tabularx}{\\columnwidth}{@{\\hskip\\tabcolsep\\extracolsep\\fill}*{' + str(int(len(table_dict)*len(GEset))) + '}{S[table-format=1.2]}}\\toprule\n') 
    for i, title in enumerate(table_dict.keys()):
        title = '\\multicolumn{2}{c}{' + title + '}'
        if 1+i<len(table_dict):
            f.write(title + ' & ')
        else:
            f.write(title + ' \\\\  \n')
    for i, title in enumerate(table_dict.keys()):
        for j, GE in enumerate(GEset):
            title = '\\multicolumn{1}{c}{' + GE + '}'
            if 1+i==len(table_dict) and 1+j == len(GEset):
                f.write(title + ' \\\\ \\midrule \n')
            else:
                f.write(title + ' & ')
    for mpcval in mpcset:
        for i, var in enumerate(table_dict.values()):
            for j, GE in enumerate(GEset):
                f.write(f"{data[mpcval]['MPCs'][GE][var]['12 months cumulative']:.2f}")
                if 1+i==len(table_dict) and 1+j == len(GEset):
                    f.write('\\\\ \n')
                else:
                    f.write(' & ')
    f.write('\\bottomrule\\end{tabularx}\n') 
    f.write('\\label{tab:' + savename + '}\n')    
    f.write('\\end{table}\n')    
    f.close()