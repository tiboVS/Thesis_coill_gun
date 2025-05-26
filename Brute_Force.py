import femm
import pandas as pd

dx = 0.1  # delta in mm
I_range = range(1, 81, 1)  # stroomsterkte iteratie

# Vasteleggen parameters projectiel en spoel
Lprojectiel = 18.7  
dProjectiel = 14.1 
DProjectiel = 16.28
StartPos = 0  

R=900

Lspoel = 90          #90         110     105   
dSpoel = 20           #20         20      16
DSpoel = 25           #25         34      24
N = 300               #300        600     300
Loop_lengte = 30     #30         20      22     lengte na de spoel tot het einde van de loop / how long the barrel is after the coil
    #130CM
    #150CM
    #90CM
M=0.1

def calc_Force(Proj_D, Proj_d, Proj_L, Proj_Start, I, x):
    femm.openfemm(1)
    femm.newdocument(0)
    femm.mi_probdef(0, 'millimeters', 'axi', 1e-8, 0, -30)  
    femm.mi_getmaterial("Air")
    femm.mi_getmaterial('Copper')  
    femm.mi_addmaterial("RVS-ferritic", 100, 100, 0, 0, 1.0e6, 0, 0, 1, 0, 0, 0)

     #boundary tekenen
    femm.mi_addboundprop('OpenBoundary', 0, 0, 1, 0, 0, 0, 0, 1)
    femm.mi_drawarc(0, -R + 100, 0, R + 100, 180, 1)
    femm.mi_drawline(0, -R + 100, 0, R + 100)
    femm.mi_setsegmentprop('OpenBoundary', 0, 1, 0, 0)
     
     
     #tekenen van spoel
    rSpoel = dSpoel/2    
    RSpoel = DSpoel/2
     
    femm.mi_drawrectangle(rSpoel,0,RSpoel,Lspoel)                              # teken een rechthoek dat de spoel voorstelt / draw a rectangle which is the coil
    femm.mi_addblocklabel((rSpoel + RSpoel) / 2, Lspoel / 2)
    femm.mi_selectlabel((rSpoel + RSpoel) / 2, Lspoel / 2)                     # plaats een lable op de spoel om het een spoel te noemen / placing lable
    femm.mi_setblockprop("Copper", 0, M, "CoilCircuit", 0, 0, N)      # wijs de materiaal parameters toe aan de spoel / tell which material it is
    femm.mi_clearselected()
    femm.mi_addblocklabel(RSpoel+2, 0)
    femm.mi_selectlabel(RSpoel+2, 0)
    femm.mi_setblockprop('Air', 1, 0, '<None>', 0, 0, 0)
    femm.mi_clearselected()
    

    # Definieer projectiel
    RProj = Proj_D / 2
    rProj = Proj_d / 2
    x1proj = Proj_Start + x
    x2proj = Proj_L + Proj_Start + x
    MiddenProj = (x1proj + x2proj) / 2

    femm.mi_drawrectangle(rProj, x1proj, RProj, x2proj)
    femm.mi_addblocklabel((rProj + RProj) / 2, MiddenProj)
    femm.mi_selectlabel((rProj + RProj) / 2, MiddenProj)
    femm.mi_setblockprop("RVS-ferritic", 0, M, "", 1, 0, 0)
    femm.mi_clearselected()

    femm.mi_addcircprop("CoilCircuit", I, 1)
    
    femm.mi_saveas("temp.fem")
    femm.mi_createmesh()
    femm.mi_analyze()
    femm.mi_loadsolution()
    femm.mo_selectblock((rProj + RProj) / 2, MiddenProj)
    F = femm.mo_blockintegral(19)
    

    femm.mi_deleteselected()
    femm.mi_clearselected()
    femm.closefemm()
    
    return F, x


for i in I_range:
    df = pd.DataFrame(columns=["Stroom", "Positie", "Kracht"])  # DataFrame per stroomsterkte
    x = -20  # Reset x per iteratie
    
    while x < 120:
        if x > 80:
            M=0.5
        B = calc_Force(DProjectiel, dProjectiel, Lprojectiel, StartPos, i, x)
        df.loc[len(df)] = [i, B[1], B[0]]
        x += dx
        

    bestandsnaam = f"data1A_i{i}.csv"
    df.to_csv(bestandsnaam, index=False, encoding="utf-8")
    print(f"Bestand opgeslagen: {bestandsnaam}")

