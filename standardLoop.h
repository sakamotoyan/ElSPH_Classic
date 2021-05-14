#pragma omp parallel for
		for (val_i i = 0; i < elFluid.numPart; i++) {
			val_i& uid = elFluid.uid(i);
			NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
			for (val_i n = 0; n < neighbs.neighbNum; n++) {
				if (isFluid(neighbs.uid[n])) {
				
				}
				else {

				}
			}
		}
#pragma omp parallel for
		for (val_i i = 0; i < elBound.numPart; i++) {
			val_i& uid = elBound.uid(i);
			NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
			for (val_i n = 0; n < neighbs.neighbNum; n++) {
				if (isFluid(neighbs.uid[n])) {
					
				}
				else {

				}
			}
		}

#pragma omp parallel for
		for (val_i i = 0; i < elFluid.numPart; i++) {
			val_i& uid = elFluid.uid(i);
			NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
			for (val_i n = 0; n < neighbs.neighbNum; n++) {
				if (isFluid(neighbs.uid[n])) {
					for (val_i ph = 0; ph < elPhase.phaseNum; ph++) {

					}
				}
			}
		}