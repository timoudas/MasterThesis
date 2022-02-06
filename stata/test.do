import delimited "C:\Users\andreas.timoudas\Python\MasterThesis\stata\data.csv", clear

TVAR_2r_grid_search sfsi bnp cpi dint, threshold(sfsi) ptrim(0.15) d(1) nlag(4) constant(1) ols(1)
