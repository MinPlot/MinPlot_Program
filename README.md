# What is MinPlot?

Information: MinPlot is a MATLAB®-based mineral formula recalculation and compositional plotting program for electron microprobe analyses (EPMA). The program offers recalculation and structural formula assignment for 15 different mineral groups: Garnet, pyroxene, olivine, amphibole, feldspar, mica, staurolite, cordierite, chlorite, chloritoid, talc, epidote, titanite, spinel, and sulfides. MinPlot is a fast and easy to use command line program and requires no prior computer programming knowledge. Percent mass fractions of oxides are loaded from datafiles and the user answers simple prompts to select mineral type, normalization scheme, and plotting options. Recalculated mineral formulas are automatically saved as output files and plots may be further manually customized by the user prior to saving. 

Loading and Saving Data: MinPlot reads data stored as text (.txt) files. The first line must contain oxide-based headers that are specific to the mineral formula to be recalculated (see Table 1 below). The headers must have capital and lowercase characters as shown in Table 1. For some phases, certain oxides are optional and will be calculated assuming a mass fraction of zero (W_j=0) if they are not included in the file read by MinPlot. MinPlot searches the header row for the column containing the appropriate header for each oxide, as a result the oxide data can be listed in any order in the input file. To start MinPlot, change the MATLAB® directory to the folder containing MinPlot and type the name of the program into the command window and hit ‘return’. When loading the data, the user is prompted to select the file in a pop-up window and, importantly, the file can be located in any folder on the user’s computer or in their network. Following calculation, the user is prompted to save their calculation. If yes, the data is automatically saved as tab delimited text files in the same directory as the source file, allowing for simplified data organization. 

<img width="758" alt="Table of Oxides" src="https://user-images.githubusercontent.com/87534196/200421489-c8585cc9-adf4-4ea5-b261-60513d7ac4c3.png">

# Mineral Plotting Options

Garnet supergroup: Garnet (X3Y2Z3O12) calculates with the following site assignments: X = Na, Ca, Ca, Mg, Mn, Fe2+, and Y at the dodecahedral site, Y = Fe3+, Cr, Ti, and viAl at the octahedral site, Z = Fe3+, ivAl, and Si at the tetrahedral site, and O2 at the anion site. Garnet structural formula are calculated using normalization to 8 cations and 12 oxygens (for Fe3+-estimation), or 12 oxygen basis alone (for totalFe=Fe2+). Endmember fractions are calculated using the matrix inversion method for solving systems of linear equations. The garnet endmembers considered are almandine (Xalm), spessartine (Xsps), grossular (Xgrs), pyrope (Xprp), andradite (Xadr), and uvarovite (Xuv). Plotting options for garnet include the Xalm + Xsps, Xgrs, and Xprp ternary. A second Fe3+, Cr, and VIAl ternary diagram for substitutions at the octahedral site is also available.

<img width="500" alt="Garnet structural formula" src="https://user-images.githubusercontent.com/87534196/200539587-0a5b4360-f02f-4920-89a0-22ae3ba148d6.png">

Pyroxene: Pyroxene (M2M1T2O6) compositions are calculated following Morimoto et al. (1989), with K, Na, Ca, Fe2+, and Mg at the distorted octahedral M2 site, Fe2+, Mg, Mn, Cr, Fe3+, Ti, and viAl at the octahedral M1 site, and Fe3+, ivAl, and Si at the tetrahedral site. Normalization is to 4 cations and 6 oxygens in the Fe3+-estimation routine, and on a 6-oxygen basis for totalFe=Fe2+. Endmember fractions are calculated using the matrix inversion method for wollastonite (Xwo), ferrosillite (Xfs), enstatite (Xen), jadeite (Xjd), aegirine (Xaeg), and kosmochlore (Xkos).

<img width="500" alt="pyroxene structural formula" src="https://user-images.githubusercontent.com/87534196/200545156-541d84d0-1006-4928-85df-e16bf07e6770.png">

A second endmember calculation option is available following Harlow (1999) where endmembers are calculated as jadeite (Xjd), aegirine (Xaeg), diopside + hedenbergite (Xdi+hd), calcium-Tschermaks pyroxene (Xcats), kosmochlor (Xkos), K-kosmochlor (XK-kos), K-jadeite (XK-jd), calcium-Eskola pyroxene (Xcaes), and enstatite (Xen).

<img width="500" alt="pyroxene structural formula" src="https://user-images.githubusercontent.com/87534196/200545427-5df23652-d60a-4342-a3bc-584b9b8f5b02.png">

Plottinng and classification follows Morimoto et al. (1989), including the 'Q-J' diagram that distinguishes Ca-Mg-Fe pyroxenes (quad), Na-Ca pyroxenes (Na-Ca), and Na pyroxenes (Na), where J = 2Na is plotted on the x-axis and Q = Ca + Mg + Fe2+ on the y-axis. For Ca-rich pyroxenes, there is a Xwo, Xfs, and Xen ternary, whereas a Xquad, Xjd, and Xaeg ternary is available for Na-rich pyroxenes. Plots of XCa-es vs Xcats and XK-cpx (XK-cpx = XK-kos + XK-jd) vs XCa-es following Harlow (1999) is available. 

Olivine: Olivine (M2TO4) is calculated here with M = Ca, Mg, Mn, Fe2+, Ni, Cr, Fe3+. Ti, and viAl at the octahedral site, and T = Fe2+, ivAl, and Si at the tetrahedral site. Normalization is to 3 cations and 4 oxygens in the Fe3+ estimation routine, and on a 4-oxygen basis for totalFe = Fe2+. Endmember fractions if forsterite (Xfo), fayalite (Xfa), tephroite (Xte), and calcio-olivine (XCa-ol) are calculated using the matrix inversion method. Plotting optiions are for the forsterite-fayalite binary, Xfo, XCa-ol, and Xfa + Xte ternary, and Fo number vs NiO (Wt %). 

<img width="400" alt="olivine structural formula" src="https://user-images.githubusercontent.com/87534196/200578690-ccc6b75c-2889-4fd0-9822-248376a75076.png">

Clinoamphibole: Amphibole (AB2C5T8O22W2) compositions are calculated following the recommendations of Leake et al. (1997) and Hawthorne et al. (2012) for structural assignment, with A = □, K, Na, and Ca at the A site, B = Ca, Na, Mn, Fe2+, and Mg at the M4 site, C = Mn, Fe2+, Mg, Fe3+, Cr, Ti, and Al at the M1, M3, and M2 sites, T = Si and Ti at the T site, and OH-, F-, Cl-, and O2- at the W site. Ferric iron is estimated following Leake et al. (1997) and Hawthorne et al. (2012), see the manuscript for details. Amphibole classification plots for Ca, Na-Ca, and Na amphiboles follow the scheme of Hawthorne et al. (2012): Amphibole compositions are plotted as C(Al + Fe3+ + 2Ti) vs A(Na + K + 2Ca). MinPlot also includes plots of Fe3+/(Al + Fe3+ + Ti) vs Fe2+/(Fe2+ + Mg + Mn) for Na amphiboles and Si vs XMg for Ca amphiboles.

Endmember fractions are not done for amphibole given the large, complicated composition space. The composition space explored for Na, Na-Ca, and Ca amphiboles in the composition diagrams is expressed by the following formulae:

<img width="800" alt="Na-amphibole structural formula" src="https://user-images.githubusercontent.com/87534196/200581309-97fa8c89-2ace-4157-81ec-01672421462e.png">

<img width="800" alt="NaCa-amphibole structural formula" src="https://user-images.githubusercontent.com/87534196/200581489-d4858ea2-4c97-42e0-95c4-efcdafb79bbc.png">

<img width="800" alt="Ca-amphibole structural formula" src="https://user-images.githubusercontent.com/87534196/200581565-f742a83e-00ba-4c95-8ccb-3686b5bbd67a.png">

Feldspar: Feldspar (AT4O8) is calculated by normalizing to 8 oxygen equivalents, with A = Ca, Na, K, Ba, Fe2+, Mn, and Mg, and T = Al and Si. Endmembers are calculated for anorthite (Xan = Ca/(Ca + Na + K)), albite (Xab = Ca/(Ca + Na + K)), and alkali feldspar (Xor = K/(Ca + Na + K)). Plotting is available as the classic An-Ab-Or feldspar ternary. An option for the subdivisions is available. 

<img width="500" alt="feldspar structural formula" src="https://user-images.githubusercontent.com/87534196/200589611-8bfe2ecb-0b52-47da-ac68-837658a2fe3b.png">

Mica: Mica (IM2-3T4O10W2) is calculated normalizing to 11 oxygen equivalents. Ions are assigned as I = □, K, Na, Ca, and Ba, M = Mg, Mn, Fe2+, Cr, Ti, and VIAl, T = IVAl and Si, and W = F, Cl, and OH. For micas totalFe is assumed to be Fe2+. MinPlot also assumes a full W site (OH = 2 – F – Cl), which may not be accurate and thus provides an estimation of the maximum possibly OH content. Mica endmembers are calculated based on two compositional groups: 1. Dioctahedral muscovite (Xms), ferroceladonite (XFe-cel), magnesioceladonite (XMg-cel), paragonite (Xpg), margarite (Xmrg), and pyrophyllite (Xprl) species, or 2. Trioctahedral, phlogopite (Xphl), annite (Xann), eastonite (Xeas), and siderophyllite (Xsid) species. The total dioctahedral or trioctahedral components are given as XDiOct and XTriOct, respectively. Plots for micas include the 1. Xms, XAl-cel, and Xprl ternary diagram, 2. Celadonite and muscovite + paragonite solid solution diagram, 3. Na (APFU) vs. Si (APFU) diagram, 4. F-Cl-OH ternary , and 5. Trioctahedral Ann-Phl-Sid-East solid solution diagram.

<img width="876" alt="mica structural formula" src="https://user-images.githubusercontent.com/87534196/200600143-7a5ae6c7-7d7f-40f3-900b-d624f34f83ce.png">


