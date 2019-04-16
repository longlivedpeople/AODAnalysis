import ROOT
import argparse



#################################################################################

def tuneHistogram(histo, color):

    histo.SetFillColorAlpha(color, 0.3)
    histo.SetLineColor(color)
    histo.GetXaxis().SetLabelSize(0.03)
    histo.GetYaxis().SetLabelSize(0.03)
    histo.GetXaxis().SetTitleSize(0.037)
    histo.GetYaxis().SetTitleSize(0.037)
    histo.GetXaxis().SetTitleOffset(1.27)
    histo.GetYaxis().SetTitleOffset(1.4)
 

#################################################################################
#################################################################################
#################################################################################


if __name__ == '__main__':

    ########## main code here
    parser = argparse.ArgumentParser(description = "Receive the parameters")
    parser.add_argument('--inputFile', action = 'store', type = str, dest = 'inputFile', help = 'Input ROOT file with the data')
    parser.add_argument('--outputDir', action = 'store', type = str, dest = 'outputDir', default = '', help = 'Name of the output dir')
    parser.add_argument('--log', action = 'store_true', dest = 'log', help = 'Linear orllogaritmic scale in Y axis')
    parser.add_argument('--norm', action = 'store', type = str, dest = 'norm', default = 'norm0', help = 'Normalization Option')
    args = parser.parse_args()


    ########## Open the tree

    file_input = ROOT.TFile(args.inputFile)

    samples = []
    variables = []

    ## List of samples and variables:
    for key in file_input.GetListOfKeys():
    
        h = key.ReadObj()
        histo_name = str(h.GetName())
        histo_name = histo_name.split('_', 1)
        if histo_name[0] not in samples: samples.append(histo_name[0])
        if histo_name[1] not in variables: variables.append(histo_name[1]) 


    ########## Load the colors
    colors = {}
    exec(open('colors.py', 'r'))


    ########## Plot the histograms

    c1 = ROOT.TCanvas('c1', '', 500, 500)
    c1.SetTicks(1,1)
    ROOT.gStyle.SetOptStat(0)
    ROOT.TGaxis.SetMaxDigits(4)
    prefix = '' # default

    for var in variables:

        ## Y axis scale options

        if (args.log): 
            prefix = 'log_'
            c1.SetLogy(1)
        else:
            prefix = 'lin_'
            c1.SetLogy(0)


        ## Legend construction

        leg = ROOT.TLegend(0.54, 0.89 - len(samples)*0.04, 0.89, 0.89)
        leg.SetTextSize(0.03)
        leg.SetBorderSize(0)


        ############# Stacked representation

        Background = ROOT.THStack("Background", "")

        for sam in samples:

            eval('file_input.' + sam + '_' + var) # Access the histogram

            if sam['type'] == 'background':

                h.SetFillColorAlpha(color[sam], 0.3)
                h.SetLineColorAlpha(ROOT.kBlack, 0.4)
                h.GetXaxis().SetLabelSize(0.03)
                h.GetYaxis().SetLabelSize(0.03)
                h.GetXaxis().SetTitleSize(0.037)
                h.GetYaxis().SetTitleSize(0.037)
                h.GetXaxis().SetTitleOffset(1.27)
                h.GetYaxis().SetTitleOffset(1.4)

                Background.Add(h)


        ############# Overlap representation implementation

        c = 0

        for sam in samples:
            
            #h = file_input.Get(sam + '_' + var)
            h = eval('file_input.' + sam + '_' + var)
            tuneHistogram(h, colors[sam])

            if (args.norm == 'norm1'):
                h.Scale(1/h.Integral())


            if (args.log): h.SetMaximum(10.*h.GetMaximum())
            else: h.SetMaximum(1.2*h.GetMaximum())

            if c == 0: 
                h.Draw('hist')
            else: 
                h.Draw('hist same')
                h.Draw('axis same')

            leg.AddEntry(h, sam, 'f')



            c+=1


        ## Draw the legend:

        leg.Draw()


        ## Save the plots

        c1.SaveAs(str(args.outputDir) + prefix + str(var) + '.png')
        c1.SaveAs(str(args.outputDir )+ prefix + str(var) + '.pdf')

        c1.Clear()


