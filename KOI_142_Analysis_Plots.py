def calc_transit_times(full_chain, specific_chains, ending_date):
    params = full_chain[specific_chains[0], specific_chains[1], :]
    filename = ("/Users/Callista/Documents/GitHub/infiles2/TTVs0" + ".in")
    infile = open(filename, 'w')
    infile.write("%.11f\n%.11f\n%.11f\n" % (g_value, M_star, params[0]))
    infile.write("%.11f %.11f %.11f %.11f %.11f %.11f\n" % (params[1], params[2], params[3], 0.0, params[4], params[5]))
    infile.write("%.11f\n" % params[6])
    infile.write("%.11f %.11f %.11f %.11f %.11f %.11f\n" % (params[7], params[8], params[9], params[11], params[10], params[12]))
    infile.close()
    
    #creating/writing setup file
    setupfilename = ("/Users/Callista/Documents/Github/setupfiles2/new_setupfile0")
    new_setupfile = open(setupfilename, 'w')
    new_setupfile.write("%s\n %.8f\n %.3f\n %d\n %d\n %d\n" % (filename, 54.675215, 0.54, ending_date, 2, 0))
    new_setupfile.close()
    os.system(("./run_TTVFast" + " " + setupfilename + " /Users/Callista/Documents/Github/KOI142_files2/final_files0" + " RV_file RV_out"))
    
    transit_times_file = np.loadtxt("/Users/Callista/Documents/Github/KOI142_files2/final_files0")
    
    planet = transit_times_file[:,0]
    epoch = transit_times_file[:,1]
    time = transit_times_file[:,2]
        
    planet_1 = planet[planet == 0]
    epoch_1 = epoch[planet == 0]
    time_1 = time[planet == 0]
    time_1 += 0./1440
    return time_1
