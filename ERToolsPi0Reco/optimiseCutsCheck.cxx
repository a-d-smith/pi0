#ifndef LARLITE_OPTIMISECUTSCHECK_CXX
#define LARLITE_OPTIMISECUTSCHECK_CXX

#include "optimiseCutsCheck.h"

namespace larlite {

  bool optimiseCutsCheck::initialize() {
		//file = new TFile("optimiseCutData.root");
		/*
		TFile file("optimiseCutData.root");
		selectedCut = (TTree*)file.Get("selectedCut");
	  cutLe = selectedCut->GetLeaf("cut");
	  evNoLe = selectedCut->GetLeaf("evNo");
	  pairNoLe = selectedCut->GetLeaf("pairNo");
		isPi0Le = selectedCut->GetLeaf("isPi0");
		file.Close();
		*/
	
		evNo = 0;
		entryNo = 0;

	  cutMin = 0.0001;
    cutMax = 0.0112;
    Ncuts  = 100; 

  	return true;
  }
  
  bool optimiseCutsCheck::analyze(storage_manager* storage) {
		auto ev_mcs = storage->get_data<event_mcshower>("mcreco"); 

		if(!ev_mcs){
      std::cout << "MCShower pointer invalid! Exiting..." << std::endl;
      exit(1);
    }

		TFile file("optimiseCutData.root");
		TTree *selectedCut = (TTree*)file.Get("selectedCut");
	  TLeaf *cutLe = selectedCut->GetLeaf("cut");
	  TLeaf *evNoLe = selectedCut->GetLeaf("evNo");
	  TLeaf *pairNoLe = selectedCut->GetLeaf("pairNo");
		TLeaf *isPi0Le = selectedCut->GetLeaf("isPi0");

    for (int n=0; n<Ncuts; n++){
      double cut = cutMin + n*(cutMax-cutMin)/Ncuts;
      int pairNo = 0;
			int i = 0;
			for (auto const &s1 : *ev_mcs){
				int j = 0;
      	for (auto const &s2 : *ev_mcs){
          if (i > j){
            if (s1.PdgCode() == 22 && s2.PdgCode() == 22){
	      			
							selectedCut->GetEntry(entryNo);
							
							// Dodgy way to make sure things match up...
							if (evNoLe->GetValue(0) == evNo){
								std::cout << cutLe->GetValue(0) - cut << std::endl;
								entryNo++;
							}
							
							pairNo++;
						}
					}
					j++;
				}
				i++;
			}
		}
    // h1->Fill(EfracMaxTolLe->GetValue(0), thetaMinTolLe->GetValue(0), case3Le->GetValue(0)/447);

		evNo++;
		file.Close();
    return true;
  }

  bool optimiseCutsCheck::finalize() {
    return true;
  }

}
#endif
