import ROOT
import sys
from array import array
import numpy as np

def calclate_Pearson_coefficient(denom_before_calib, num_PKU, denom_after_calib, var_type):
    light_yield_denom_mean = np.mean(denom_before_calib)
    light_yield_num_mean = np.mean(num_PKU)
    light_yield_denom_calib_mean = np.mean(denom_after_calib)

    Pearson_denom_1 = (
        np.sqrt(sum([(x-light_yield_denom_mean)**2/(len(denom_before_calib)-1) for x in denom_before_calib])) * 
        np.sqrt(sum([(x-light_yield_num_mean)**2/(len(num_PKU)-1) for x in num_PKU]))
        )
    Pearson_num_1 = (
        sum([(x-light_yield_denom_mean)*(y-light_yield_num_mean)/(len(denom_before_calib)-1) for x,y in zip(denom_before_calib, num_PKU) ])
    )

    print(f"Before calib, {var_type} Pearson coefficient = ", Pearson_num_1/Pearson_denom_1, "----------------------------------")

    Pearson_denom_2 = (
        np.sqrt(sum([(x-light_yield_denom_calib_mean)**2/(len(denom_before_calib)-1) for x in denom_after_calib])) * 
        np.sqrt(sum([(x-light_yield_num_mean)**2/(len(num_PKU)-1) for x in num_PKU]))
        )
    Pearson_num_2 = (
        sum([(x-light_yield_denom_calib_mean)*(y-light_yield_num_mean)/(len(denom_before_calib)-1) for x,y in zip(denom_after_calib, num_PKU) ])
    )
    print(f"After calib, {var_type} Pearson coefficient = ", Pearson_num_2/Pearson_denom_2, "----------------------------------")


def calib_for_one(file1_name, file2_name, output_file_name):
    file1 = ROOT.TFile.Open(file1_name)
    file2 = ROOT.TFile.Open(file2_name)

    output_file = ROOT.TFile(output_file_name, "RECREATE")

    graph_list = ["g_spe_L_vs_bar", "g_spe_R_vs_bar", "g_lyso_L_pc_per_kev_vs_bar", "g_lyso_R_pc_per_kev_vs_bar", "g_L_light_yield_vs_bar", "g_R_light_yield_vs_bar"]
    light_yield_g = ["g_L_light_yield_vs_bar", "g_R_light_yield_vs_bar"]
    ratio_g_list = []
    
    spe_lyso_L_R_denom_num = [[],[],[],[],[],[],[],[],[],[],[],[]]

    for _graph in graph_list:


        graph1 = file1.Get(f"{_graph}")
        graph2 = file2.Get(f"{_graph}")

        if not graph1 or not graph2:
            print("Error: cannot read from th Graph")
            exit(1)

    
        n_points1 = graph1.GetN()
        n_points2 = graph2.GetN()

        if n_points1 != n_points2:
            print("Error: have different dot number")
            exit(1)

        # create new graph
        ratio_graph = ROOT.TGraph()
        ratio_graph.SetName(f"{_graph}_ratio") #name will be used in another function
        
        if "spe" in _graph or "lyso" in _graph:
            ratio_g_list.append(f"{_graph}_ratio")

        # run all dots
        for i in range(n_points1):
            x1 = array('d', [0])
            x2 = array('d', [0])
            y1 = array('d', [0])
            y2 = array('d', [0])
        
            graph1.GetPoint(i, x1, y1)
            graph2.GetPoint(i, x2, y2)

            if abs(x1[0] - x2[0]) > 1e-6:
                print(f"Warning: value of x of {i}th dot is different: {x1} vs {x2}, jump this dot")
                continue

            if y2[0] != 0:
                ratio = y1[0] / y2[0]
            else:
                ratio = 0 # set NAN or float('nan')

            ratio_graph.SetPoint(i, x1[0], ratio)

            if "spe_L" in _graph:
                spe_lyso_L_R_denom_num[0].append(y2[0])
                spe_lyso_L_R_denom_num[1].append(y1[0])

            if "spe_R" in _graph:
                spe_lyso_L_R_denom_num[2].append(y2[0])
                spe_lyso_L_R_denom_num[3].append(y1[0])

            if "lyso_L" in _graph:
                spe_lyso_L_R_denom_num[4].append(y2[0])
                spe_lyso_L_R_denom_num[5].append(y1[0])

            if "lyso_R" in _graph:
                spe_lyso_L_R_denom_num[6].append(y2[0])
                spe_lyso_L_R_denom_num[7].append(y1[0])

            if "L_light" in _graph:
                spe_lyso_L_R_denom_num[8].append(y2[0])
                spe_lyso_L_R_denom_num[9].append(y1[0])

            if "R_light" in _graph:
                spe_lyso_L_R_denom_num[10].append(y2[0])
                spe_lyso_L_R_denom_num[11].append(y1[0])
        
        ratio_graph.Write()

    # save to the new root file
    output_file.Close()

    print(f"Completed: have saved graph_ratio to {output_file_name}")

    return ratio_g_list, spe_lyso_L_R_denom_num, light_yield_g


def average(file_list, graph_list, output_file_name, spe_lyso_L_R_denom_num, light_yield_g):

    Rfile = []
    
    # create a file list
    for _file in file_list:
        Rfile.append(ROOT.TFile.Open(_file))

    output_file = ROOT.TFile(output_file_name, "RECREATE")

    #save the spe average sf in the list for ratio for all channels
    for_ratio_spe_L = []
    for_ratio_spe_R = []
    #same for lyso
    for_ratio_lyso_L = []
    for_ratio_lyso_R = []

    #save the spe not for all just for average
    for_ratio_spe_L_not_all = []
    spe_L_std_dev = []
    for_ratio_spe_R_not_all = []
    spe_R_std_dev = []
    #same for lyso
    for_ratio_lyso_L_not_all = []
    lyso_L_std_dev = []
    for_ratio_lyso_R_not_all = []
    lyso_R_std_dev = []

    for _graph in graph_list:
        
        Rgraph = []
        #create a TGraph list
        for _Rfile in Rfile:
            Rgraph.append(_Rfile.Get(f'{_graph}'))
        
        points = Rgraph[0].GetN()
        
        #average graph
        ave_graph = ROOT.TGraphErrors()
        ave_graph.SetName(f'{_graph}_ave')
        ave_list = []
        ave_std_dev_list = []
        ave_list_all = []

        for_spe_lyso = ROOT.TGraph()
        for_spe_lyso_calib = ROOT.TGraph()

        fit_func_before = ROOT.TF1("fit_func_before", "[0]*x", 3.4, 4.4)
        fit_func_after = ROOT.TF1("fit_func_after", "[0]*x", 3.4, 4.4)
        
        # loop all dots (default 16 channels)
        for i in range(points):
            x = []
            y = []

            # reset x and y for each loop
            for j in range(len(Rgraph)):
                x.append(array('d', [0]))
                y.append(array('d', [0]))
            
            #get point for each TGraph
            for j_R in range(len(Rgraph)):
                Rgraph[j_R].GetPoint(i, x[j_R], y[j_R])

            ave = 0
            sq_dev = 0
            std_dev = 0
            
            #according to the Nb of TGraph calculate average value
            for j_a in range(len(Rgraph)):
                ave = ave + y[j_a][0]/len(Rgraph)
            
            #calculate error bar
            for j_a_1 in range(len(Rgraph)):
                sq_dev = sq_dev + (ave-y[j_a_1][0])**2/(len(Rgraph)-1)
            
            std_dev = np.sqrt(sq_dev)
            
            ave_list.append(ave)
            ave_std_dev_list.append(std_dev)
            ave_graph.SetPoint(i, x[0][0], ave)
            ave_graph.SetPointError(i, 0, std_dev)

        #make average list length equal to all channels, 
        #and each element is the calibartion value for corresponding channel 
        for j_a_2 in range(len(Rgraph)):
            ave_list_all.extend(ave_list)

        #scatter plots for before calib and after calib
        if "spe_L" in _graph:
            for_ratio_spe_L = ave_list_all
            for_ratio_spe_L_not_all = ave_list
            spe_L_std_dev = ave_std_dev_list

            spe_L_denom = array('d', spe_lyso_L_R_denom_num[0])
            spe_L_denom_calib = array('d',[x*y for x, y in zip(spe_L_denom, ave_list_all)])
            spe_L_num = array('d', spe_lyso_L_R_denom_num[1])
            
            for_spe_lyso = ROOT.TGraph(len(spe_lyso_L_R_denom_num[0]), spe_L_denom, spe_L_num)
            for_spe_lyso.SetName("spe_L")
            for_spe_lyso_calib = ROOT.TGraph(len(spe_lyso_L_R_denom_num[0]), spe_L_denom_calib, spe_L_num)
            for_spe_lyso_calib.SetName("spe_L_calib")
            
            print("Before calib, dev for spe_L is = ", np.sqrt(sum([(x-y)**2 for x,y in zip(spe_L_denom, spe_L_num)])))
            print("After calib, dev for spe_L is = ", np.sqrt(sum([(x-y)**2 for x,y in zip(spe_L_denom_calib, spe_L_num)])))

            calclate_Pearson_coefficient(spe_L_denom, spe_L_num, spe_L_denom_calib, "spe_L")

            for_spe_lyso.Fit(fit_func_before)
            fit_func_before.SetName("spe_L_linear_fit_func_before_calib")
            for_spe_lyso_calib.Fit(fit_func_after)
            fit_func_after.SetName("spe_L_linear_fit_func_after_calib")

        if "spe_R" in _graph:
            for_ratio_spe_R = ave_list_all
            for_ratio_spe_R_not_all = ave_list
            spe_R_std_dev = ave_std_dev_list

            spe_R_denom = array('d', spe_lyso_L_R_denom_num[2])
            spe_R_num = array('d', spe_lyso_L_R_denom_num[3])
            spe_R_denom_calib = array('d',[x*y for x, y in zip(spe_R_denom, ave_list_all)])
            
            for_spe_lyso = ROOT.TGraph(len(spe_lyso_L_R_denom_num[2]), spe_R_denom, spe_R_num)
            for_spe_lyso.SetName("spe_R")
            for_spe_lyso_calib = ROOT.TGraph(len(spe_lyso_L_R_denom_num[2]), spe_R_denom_calib, spe_R_num)
            for_spe_lyso_calib.SetName("spe_R_calib")
            
            print("Before calib, dev for spe_R is = ", np.sqrt(sum([(x-y)**2 for x,y in zip(spe_R_denom, spe_R_num)])))
            print("After calib, dev for spe_R is = ", np.sqrt(sum([(x-y)**2 for x,y in zip(spe_R_denom_calib, spe_R_num)])))

            calclate_Pearson_coefficient(spe_R_denom, spe_R_num, spe_R_denom_calib, "spe_R")

            for_spe_lyso.Fit(fit_func_before)
            fit_func_before.SetName("spe_R_linear_fit_func_before_calib")
            for_spe_lyso_calib.Fit(fit_func_after)
            fit_func_after.SetName("spe_R_linear_fit_func_after_calib")
   
        if "lyso_L" in _graph:
            for_ratio_lyso_L = ave_list_all
            for_ratio_lyso_L_not_all = ave_list
            lyso_L_std_dev = ave_std_dev_list

            lyso_L_denom = array('d', spe_lyso_L_R_denom_num[4])
            lyso_L_num = array('d', spe_lyso_L_R_denom_num[5])
            lyso_L_denom_calib = array('d',[x*y for x, y in zip(lyso_L_denom, ave_list_all)])
            
            for_spe_lyso = ROOT.TGraph(len(spe_lyso_L_R_denom_num[4]), lyso_L_denom, lyso_L_num)
            for_spe_lyso.SetName("lyso_L")
            for_spe_lyso_calib = ROOT.TGraph(len(spe_lyso_L_R_denom_num[4]), lyso_L_denom_calib, lyso_L_num)
            for_spe_lyso_calib.SetName("lyso_L_calib")
            
            print("Before calib, dev for lyso_L is = ", np.sqrt(sum([(x-y)**2 for x,y in zip(lyso_L_denom, lyso_L_num)])))
            print("After calib, dev for lyso_L is = ", np.sqrt(sum([(x-y)**2 for x,y in zip(lyso_L_denom_calib, lyso_L_num)])))

            calclate_Pearson_coefficient(lyso_L_denom, lyso_L_num, lyso_L_denom_calib, "lyso_L")

            for_spe_lyso.Fit(fit_func_before)
            fit_func_before.SetName("lyso_L_linear_fit_func_before_calib")
            for_spe_lyso_calib.Fit(fit_func_after)
            fit_func_after.SetName("lyso_L_linear_fit_func_after_calib")

        if "lyso_R" in _graph:
            for_ratio_lyso_R = ave_list_all
            for_ratio_lyso_R_not_all = ave_list
            lyso_R_std_dev = ave_std_dev_list

            lyso_R_denom = array('d', spe_lyso_L_R_denom_num[6])
            lyso_R_num = array('d', spe_lyso_L_R_denom_num[7])
            lyso_R_denom_calib = array('d',[x*y for x, y in zip(lyso_R_denom, ave_list_all)])
            
            for_spe_lyso = ROOT.TGraph(len(spe_lyso_L_R_denom_num[6]), lyso_R_denom, lyso_R_num)
            for_spe_lyso.SetName("lyso_R")
            for_spe_lyso_calib = ROOT.TGraph(len(spe_lyso_L_R_denom_num[6]), lyso_R_denom_calib, lyso_R_num)
            for_spe_lyso_calib.SetName("lyso_R_calib")
            
            print("Before calib, dev for lyso_R is = ", np.sqrt(sum([(x-y)**2 for x,y in zip(lyso_R_denom, lyso_R_num)])))
            print("After calib, dev for lyso_R is = ", np.sqrt(sum([(x-y)**2 for x,y in zip(lyso_R_denom_calib, lyso_R_num)])))

            calclate_Pearson_coefficient(lyso_R_denom, lyso_R_num, lyso_R_denom_calib, "lyso_R")

            for_spe_lyso.Fit(fit_func_before)
            fit_func_before.SetName("lyso_R_linear_fit_func_before_calib")
            for_spe_lyso_calib.Fit(fit_func_after)
            fit_func_after.SetName("lyso_R_linear_fit_func_after_calib")
                
 
        ave_graph.Write()
        for_spe_lyso.Write()
        for_spe_lyso_calib.Write()
        fit_func_before.Write()
        fit_func_after.Write()
    
    for _graph in light_yield_g:
        for_light_yield = ROOT.TGraph()
        for_light_yield_calib = ROOT.TGraph()
        ave_for_LY = ROOT.TGraphErrors()

        fit_func_before_ly = ROOT.TF1("fit_func_before", "[0]*x", 2400, 3500)
        fit_func_after_ly = ROOT.TF1("fit_func_after", "[0]*x", 2400, 3500)
        if "L" in _graph:
            spe_lyso_ratio = array('d', [x/y for x,y in zip(for_ratio_lyso_L, for_ratio_spe_L)])
            light_yield_denom = array('d', spe_lyso_L_R_denom_num[8])
            light_yield_num = array('d', spe_lyso_L_R_denom_num[9])
            light_yield_denom_calib = array('d', [x*y for x,y in zip(light_yield_denom, spe_lyso_ratio)])
                
            for_light_yield = ROOT.TGraph(len(spe_lyso_L_R_denom_num[8]), light_yield_denom, light_yield_num)
            for_light_yield.SetName("light_yield_L")
            for_light_yield_calib = ROOT.TGraph(len(spe_lyso_L_R_denom_num[8]), light_yield_denom_calib, light_yield_num)
            for_light_yield_calib.SetName("light_yield_L_calib")

            print("Before calib, dev for light yield left is = ", np.sqrt(sum([(x-y)**2 for x,y in zip(light_yield_denom, light_yield_num)])))
            print("After calib, dev for light yield left is = ", np.sqrt(sum([(x-y)**2 for x,y in zip(light_yield_denom_calib, light_yield_num)])))

            calclate_Pearson_coefficient(light_yield_denom, light_yield_num, light_yield_denom_calib, "light_yield_L")

            for_light_yield.Fit(fit_func_before_ly)
            fit_func_before_ly.SetName("light_yield_L_before_calib")
            for_light_yield_calib.Fit(fit_func_after_ly)
            fit_func_after_ly.SetName("light_yield_L_after_calib")

            for i in range(len(for_ratio_spe_L_not_all)):
                pass_error = np.sqrt((spe_L_std_dev[i]/for_ratio_spe_L_not_all[i])**2 + (lyso_L_std_dev[i]/for_ratio_lyso_L_not_all[i]))
                ratio = for_ratio_lyso_L_not_all[i]/for_ratio_spe_L_not_all[i]
                ave_for_LY.SetPoint(i, i, ratio)
                ave_for_LY.SetPointError(i, 0, pass_error)
                ave_for_LY.SetName("SF_for_light_yield_L")
        
        if "R" in _graph:
            spe_lyso_ratio = array('d', [x/y for x,y in zip(for_ratio_lyso_R, for_ratio_spe_R)])
            light_yield_denom = array('d', spe_lyso_L_R_denom_num[10])
            light_yield_num = array('d', spe_lyso_L_R_denom_num[11])
            light_yield_denom_calib = array('d', [x*y for x,y in zip(light_yield_denom, spe_lyso_ratio)])
                
            for_light_yield = ROOT.TGraph(len(spe_lyso_L_R_denom_num[10]), light_yield_denom, light_yield_num)
            for_light_yield.SetName("light_yield_R")
            for_light_yield_calib = ROOT.TGraph(len(spe_lyso_L_R_denom_num[10]), light_yield_denom_calib, light_yield_num)
            for_light_yield_calib.SetName("light_yield_R_calib")

            print("Before calib, dev for light yield right is = ", np.sqrt(sum([(x-y)**2 for x,y in zip(light_yield_denom, light_yield_num)])))
            print("After calib, dev for light yield right is = ", np.sqrt(sum([(x-y)**2 for x,y in zip(light_yield_denom_calib, light_yield_num)])))

            calclate_Pearson_coefficient(light_yield_denom, light_yield_num, light_yield_denom_calib, "light_yield_R")

            for_light_yield.Fit(fit_func_before_ly)
            fit_func_before_ly.SetName("light_yield_R_before_calib")
            for_light_yield_calib.Fit(fit_func_after_ly)
            fit_func_after_ly.SetName("light_yield_R_after_calib")

            for i in range(len(for_ratio_spe_R_not_all)):
                pass_error = np.sqrt((spe_R_std_dev[i]/for_ratio_spe_R_not_all[i])**2 + (lyso_R_std_dev[i]/for_ratio_lyso_R_not_all[i]))
                ratio = for_ratio_lyso_R_not_all[i]/for_ratio_spe_R_not_all[i]
                ave_for_LY.SetPoint(i, i, ratio)
                ave_for_LY.SetPointError(i, 0, pass_error)
                ave_for_LY.SetName("SF_for_light_yield_R")
        
        ave_for_LY.Write()
        for_light_yield.Write()
        for_light_yield_calib.Write()
        fit_func_before_ly.Write()
        fit_func_after_ly.Write()

    
    output_file.Close()

    print("calculate average Done !!!")

def for_each_BAC(B_list, B_name):
    
    spe_lyso_L_R_denom_num = [[],[],[],[],[],[],[],[],[],[],[],[]]

    file_list = []
    graph_list = []
    light_yield_list = []
    for _B in range(len(B_list)):
        _spe_lyso_L_R_denom_num = [[],[],[],[],[],[],[],[],[],[],[],[]]
        graph_list, _spe_lyso_L_R_denom_num, light_yield_list= calib_for_one(f"../{B_name}/in_PKU/module_3211002000{B_list[_B]}_analysis.root\
                      ", f"../{B_name}/in_{B_name}/module_3211002000{B_list[_B]}_analysis.root\
                        ", f"../{B_name}/out_3211002000{B_list[_B]}.root")
        file_list.append(f"../{B_name}/out_3211002000{B_list[_B]}.root")

        for i in range(len(spe_lyso_L_R_denom_num)):
            spe_lyso_L_R_denom_num[i].extend(_spe_lyso_L_R_denom_num[i])

    average(file_list, graph_list, f"ave_{B_name}.root", spe_lyso_L_R_denom_num, light_yield_list)
    print(f"Calibration for {B_name} is DONE!!! --------------- file saved as ./ave_{B_name}.root, please check!!!")

def main():
    MIB = ['0518', '0545', '0531', '0547']# barcodes of SMs test at MIB
    UVA = ['5928', '5952','0518', '0545', '0531', '0547'] # barcodes of SMs test at UVA
    CIT = ['8421', '8419', '8418', '8621', '8519'] # barcodes of SMs test at CIT
    PKU = ['2826', '2952', '2934', '3097'] # barcodes of SMs being retested in PKU
    
    # calculating average value for Milano
    for_each_BAC(MIB, "MIB")
    for_each_BAC(UVA, "UVA")
    for_each_BAC(CIT, "CIT")
    
    # for PKU to do calib for different periods
    # spe_lyso_L_R_denom_num = [[],[],[],[],[],[],[],[],[],[],[],[]]

    # file_list = []
    # graph_list = []
    # light_yield_list = []
    # for _B in range(len(PKU)):
    #     _spe_lyso_L_R_denom_num = [[],[],[],[],[],[],[],[],[],[],[],[]]
    #     graph_list, _spe_lyso_L_R_denom_num, light_yield_list= calib_for_one(f"../PKU/new/module_3211002000{PKU[_B]}_analysis.root\
    #                   ", f"../PKU/old/module_3211002000{PKU[_B]}_analysis.root\
    #                     ", f"../PKU/out_3211002000{PKU[_B]}.root")
    #     file_list.append(f"../PKU/out_3211002000{PKU[_B]}.root")

    #     for i in range(len(spe_lyso_L_R_denom_num)):
    #         spe_lyso_L_R_denom_num[i].extend(_spe_lyso_L_R_denom_num[i])

    # average(file_list, graph_list, f"ave_PKU.root", spe_lyso_L_R_denom_num, light_yield_list)

main()