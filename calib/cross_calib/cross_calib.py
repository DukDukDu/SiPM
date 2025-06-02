import ROOT
import sys
from array import array

def calib_for_one(file1_name, file2_name, output_file_name):
    file1 = ROOT.TFile.Open(file1_name)
    file2 = ROOT.TFile.Open(file2_name)

    output_file = ROOT.TFile(output_file_name, "RECREATE")

    graph_list = ["g_spe_L_vs_bar", "g_spe_R_vs_bar", "g_lyso_L_pc_per_kev_vs_bar", "g_lyso_R_pc_per_kev_vs_bar"]
    ratio_g_list = []

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
        
        ratio_graph.Write()

    # save to the new root file
    output_file.Close()

    print(f"Completed: have saved graph_ratio to {output_file_name}")

    return ratio_g_list


def average(file_list, graph_list, output_file_name):

    Rfile = []
    
    # create a file list
    for _file in file_list:
        Rfile.append(ROOT.TFile.Open(_file))

    output_file = ROOT.TFile(output_file_name, "RECREATE")

    for _graph in graph_list:
        Rgraph = []
        #create a TGraph list
        for _Rfile in Rfile:
            Rgraph.append(_Rfile.Get(f'{_graph}'))
        
        points = Rgraph[0].GetN()
        
        #average graph
        ave_graph = ROOT.TGraph()
        ave_graph.SetName(f'{_graph}_ave')
        
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
            
            #according to the Nb of TGraph calculate average value
            for j_a in range(len(Rgraph)):
                ave = ave + y[j_a][0]/len(Rgraph)

            ave_graph.SetPoint(i, x[0][0], ave)

        ave_graph.Write()

    output_file.Close()

    print("calculate average Done !!!")

def for_each_BAC(B_list, B_name):
    file_list = []
    graph_list = []
    for _B in range(len(B_list)):
        graph_list = calib_for_one(f"../{B_name}/in_PKU/module_3211002000{B_list[_B]}_analysis.root\
                      ", f"../{B_name}/in_{B_name}/module_3211002000{B_list[_B]}_analysis.root\
                        ", f"../{B_name}/out_3211002000{B_list[_B]}.root")
        file_list.append(f"../{B_name}/out_3211002000{B_list[_B]}.root")

    average(file_list, graph_list, f"ave_{B_name}.root")
    print(f"Calibration for {B_name} is DONE!!! --------------- file saved as ./ave_{B_name}.root, please check!!!")

def main():
    MIB = []# barcodes of SMs test at MIB
    UVA = ['5928', '5952','0518', '0545', '0531', '0547'] # barcodes of SMs test at UVA
    CIT = ['8421', '8419', '8418', '8621', '8519'] # barcodes of SMs test at CIT
    PKU_retest = ['2826', '2952', '2934', '3097'] # barcodes of SMs being retested in PKU
    
    # calculating average value for Milano
    for_each_BAC(MIB, "MIB")
    for_each_BAC(UVA, "UVA")
    for_each_BAC(CIT, "CIT")
    # for_each_BAC(CIT, "CIT")
    
    
    
    # for PKU to do calib for different period
    file_list = []
    graph_list = []
    for _P in range(len(PKU_retest)):
        graph_list = calib_for_one(f"../PKU/new/module_3211002000{PKU_retest[_P]}_analysis.root\
                      ", f"../PKU/old/module_3211002000{PKU_retest[_P]}_analysis.root\
                        ", f"../PKU/out_3211002000{PKU_retest[_P]}.root")
        file_list.append(f"../PKU/out_3211002000{PKU_retest[_P]}.root")

    average(file_list, graph_list, "ave_PKU.root")

main()