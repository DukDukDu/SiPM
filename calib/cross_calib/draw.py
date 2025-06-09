import ROOT

file = ROOT.TFile.Open("ave_MIB.root")

graph = file.Get("light_yield_R_calib")
if not graph:
    raise RuntimeError("cannot find TGraphErrors ")

canvas = ROOT.TCanvas("canvas", "TGraphErrors Example", 800, 600)


# graph.GetXaxis().SetLimits(-1, 16)    
# graph.SetMinimum(0)                    
# graph.SetMaximum(1.2)                  

graph.SetTitle(";LY_R tested in PKU;LY_R tested in MIB")
graph.SetMarkerStyle(20)
graph.SetMarkerColor(ROOT.kBlue)


graph.Draw("AP")  # A: axis, P: point (with error bars)


# canvas.Draw()


canvas.SaveAs("SF_graph.png")
