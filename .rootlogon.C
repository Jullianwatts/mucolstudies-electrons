{
    printf("Running rootlogon.C\n");
    //gROOT->LoadMacro("/scratch/trholmes/utils/tdrstyle.C"); 
    //setTDRStyle();
    gROOT->LoadMacro("/scratch/jwatts/mucol/utils/AtlasStyle.C");
    SetAtlasStyle();
}
