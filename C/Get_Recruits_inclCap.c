void Get_Recruits(MSEBoxModel *bm, int species, int stock_id, double plankton, FILE *llogfp) {
	double step1, step2, jack_SSB, jack_B, jack_a, larval_scalar;
	double temprec = 0.0;
	double testrec = 0.0;
	double testrec20 = 0.0;
	double testrech = 0.0;
	double testrecap = 0.0;
	double BHalpha_sp = FunctGroupArray[species].speciesParams[BHalpha_id];
	double BHbeta_sp = FunctGroupArray[species].speciesParams[BHbeta_id];
	double Ralpha_sp = FunctGroupArray[species].speciesParams[Ralpha_id];
	double Rbeta_sp = FunctGroupArray[species].speciesParams[Rbeta_id];
	double recover_mult_sp = FunctGroupArray[species].speciesParams[recover_mult_id];
	double recover_start_sp = FunctGroupArray[species].speciesParams[recover_start_id];
	//double KWSR_sp = FunctGroupArray[species].speciesParams[KWSR_id];
	//double KWRR_sp = FunctGroupArray[species].speciesParams[KWRR_id];
	int recruit_sp = (int) (FunctGroupArray[species].speciesParams[flagrecruit_id]);
	int sp_numGeneTypes = (int) (FunctGroupArray[species].numGeneTypes);
	int ngene = 0;
	int do_debug = 0;
    int qid = EMBRYO[species].next_spawn;

    if(EMBRYO[species].SpawnRecruitOverlap)  // Do not use the EMBRYO[species].CounterNotDone chekck here as it throws things out of whack for some species who spawn during teh same period recruits arrive
        qid++;
    
    if(((bm->current_box == bm->checkbox) || (bm->checkbox > bm->nbox)) && (bm->dayt >= bm->checkstart) && (bm->which_check == species))
		do_debug = 1;
    
    if((bm->dayt >= bm->checkstart) && (bm->which_check == species))
		do_debug = 1;
    
    if ((FunctGroupArray[species].numSpawns > 1) || (FunctGroupArray[species].groupType == FISH_INVERT))
        do_debug = 1;

    //do_debug = 1;
    
	for(ngene = 0; ngene < sp_numGeneTypes; ngene++){
		switch (recruit_sp) {
		case no_recruit:
			quit("No such flagrecruit defined (i.e. value must be > 0)\n");
			break;
		case const_recruit: /* Fixed set of constants */
		case chl_recruit: /* Proportional to primary productivity */
		case rand_recruit: /* Random - follows lognormal */
		case plank_recruit: /* Spawn is based on plankton levels (not just CHLa) */
		case linear_recruit:/* Pupping or calving linearly dependent on maternal condition */
		case fixed_linear_recruit:/* Pupping or calving a fixed number per adult spawning */
		case ts_recruit:/* Read in timeseries of recruitment */
        case multiple_ts_recruit: /* Read in timeseries of recruitment */
			temprec = EMBRYO[species].Larvae[stock_id][ngene][qid];
			break;
		case BevHolt_recruit: /* Beverton-Holt stock-recruit relationship - Atlantis basic version (mix numbers and biomass) */
			temprec = (recSTOCK[species][stock_id] * BHalpha_sp * EMBRYO[species].Larvae[stock_id][ngene][qid] / (BHbeta_sp + bm->totfishpop[species]
					* stock_prop[species][stock_id]));

            /**/
			if (do_debug && (bm->which_check == species)) {
				fprintf( llogfp,
					"Time: %e, box%d-%d, species %s, ngene: %d, stock: %d, temprec: %e, recSTOCK: %e, BHalpha: %e, Larvae: %e, BHbeta: %e, totfish: %e, stock_prop: %e)\n",
						bm->dayt, bm->current_box, bm->current_layer, FunctGroupArray[species].groupCode, 
						ngene, stock_id, temprec, recSTOCK[species][stock_id], BHalpha_sp,
						EMBRYO[species].Larvae[stock_id][ngene][qid], BHbeta_sp, bm->totfishpop[species],
						stock_prop[species][stock_id]);
			}
             /**/
			break;
        case BevHolt_num_recruit: /* Beverton-Holt stock-recruit relationship - numbers only case */
            temprec = (recSTOCK[species][stock_id] * BHalpha_sp * EMBRYO[species].Larvae[stock_id][ngene][qid] / (BHbeta_sp  + recSTOCK[species][stock_id] * EMBRYO[species].Larvae[stock_id][ngene][qid]));

            break;
		case BevHolt_rand_recruit: /* Spawn is based on Beverton Holt with lognormal variation and
		 dependence on plankton levels */
			step1 = Util_Logx_Result(-lognorm_mu, lognorm_sigma);
			step2 = (recSTOCK[species][stock_id] * BHalpha_sp * EMBRYO[species].Larvae[stock_id][ngene][qid] / (BHbeta_sp + bm->totfishpop[species]
					* stock_prop[species][stock_id])) * (plankton / bm->ref_chl);
			temprec = step2 * step1;
			break;
		case recover_recruit: /* Spawn is allowed a recovery encouraging boost of recruits
		 after "recovery_span" years of depressed stock levels */
			temprec = (recSTOCK[species][stock_id] * BHalpha_sp * EMBRYO[species].Larvae[stock_id][ngene][qid] / (BHbeta_sp + bm->totfishpop[species]
					* stock_prop[species][stock_id]));

			if (recover_help[species][0] && (recover_help_set[species] <= bm->dayt - (recover_span * recover_subseq))) {
				temprec *= recover_mult_sp;
				if (recover_help[species][1]) {
					fprintf(llogfp, "Time: %e, species %s has had a recovery event (temprec = %e)\n", bm->dayt, FunctGroupArray[species].groupCode, temprec);
					recover_help[species][1] = 0;
				}
			}
			break;
		case force_recover_recruit: /* Spawn has a pre-specified recovery encouraging boost of recruits */
			temprec = (recSTOCK[species][stock_id] * BHalpha_sp * EMBRYO[species].Larvae[stock_id][ngene][qid] / (BHbeta_sp + bm->totfishpop[species]
					* stock_prop[species][stock_id]));

			if ((bm->dayt >= recover_start_sp) && (bm->dayt <= (recover_start_sp + recover_subseq))) {
				temprec *= recover_mult_sp;
				fprintf(llogfp, "Time: %e, species %s has had a prescribed recovery event (temprec = %e)\n", bm->dayt, FunctGroupArray[species].groupCode, temprec);
			}
			break;
		case Ricker_recruit: /* Ricker */
			temprec = bm->totfishpop[species] * stock_prop[species][stock_id] * exp(recSTOCK[species][stock_id] * Ralpha_sp * (1.0 - bm->totfishpop[species] * stock_prop[species][stock_id] / Rbeta_sp));
            /*
            if (do_debug && (bm->which_check == species)) {
                fprintf( llogfp, "Time: %e, box%d-%d, species %s, ngene: %d, stock: %d, temprec: %e, totfishpop: %e, stock_prop: %e, recSTOCK: %e, Ralpha: %e, Larvae: %e, Rbeta: %e\n",
                    bm->dayt, bm->current_box, bm->current_layer, FunctGroupArray[species].groupCode,
                    ngene, stock_id, temprec, bm->totfishpop[species], stock_prop[species][stock_id], recSTOCK[species][stock_id], Ralpha_sp,
                    EMBRYO[species].Larvae[stock_id][ngene][qid], Rbeta_sp);
            }
             */
                
			break;
		case SSB_BevHolt_recruit: /* Standard Beverton-Holt stock-recruit relationship */
			temprec = (recSTOCK[species][stock_id] * BHalpha_sp * bm->totfishpop[species] * stock_prop[species][stock_id] / (BHbeta_sp + bm->totfishpop[species] * stock_prop[species][stock_id]));

			testrec = ((1 / 0.9) * recSTOCK[species][stock_id] * stock_prop[species][stock_id] * BHalpha_sp * bm->totinitpop[species]) / (BHbeta_sp + (1 / 0.9) *  stock_prop[species][stock_id] * bm->totinitpop[species]);
			testrec20 = ((0.2 / 0.9) * recSTOCK[species][stock_id] * stock_prop[species][stock_id] * BHalpha_sp * bm->totinitpop[species]) / (BHbeta_sp + (0.2 / 0.9) *  stock_prop[species][stock_id] * bm->totinitpop[species]);
			testrech = testrec20 / testrec; //calculated h - initial conditions at 0.9 B0, hence the 0.9 above
			testrecap = testrec * (1 + (testrech / 50)); //50 too large here
			/*if (temprec > testrecap) { // alt using fn of R0
				temprec = testrecap; //caps at a fn of R0 and h
			}*/
			if (temprec > testrec) {
				temprec = testrec; //caps at R0 instead of cap which is a fn of R0 and h
			}
			/*if (bm->which_check == species) {
				fprintf(llogfp, "Time: %e, species %s, testrec %e, testrec20 %e, testrech %e, restrecap %e, alpha %e, temprec %e, thisBiomass %e, beta %e\n", bm->dayt, FunctGroupArray[species].groupCode, testrec, testrec20, testrech, testrecap, BHalpha_sp, temprec, bm->totfishpop[species], BHbeta_sp);
			}*/
			break;
		case jackknife_recruit: /* Jackknife spawning function */
			jack_SSB = EMBRYO[species].Larvae[stock_id][ngene][qid];
			jack_B = FunctGroupArray[species].speciesParams[jack_b_id] * bm->totfishpop[species] * stock_prop[species][stock_id];
			jack_a = FunctGroupArray[species].speciesParams[jack_a_id];

			if (jack_SSB <  jack_B)
				temprec = jack_a * jack_SSB;
			else
				temprec = jack_a * jack_B;
			break;
		case baltic_ricker: /* Baltic version of the ricker */
            temprec = 1000 * recSTOCK[species][stock_id] * Ralpha_sp * bm->totfishpop[species] * bm->X_CN * mg_2_tonne
                * stock_prop[species][stock_id] * exp( -1.0 * Rbeta_sp * bm->totfishpop[species] * bm->X_CN * mg_2_tonne * stock_prop[species][stock_id]);
			break;
        case SSB_ricker: /* SSB based version of the ricker given short lived and senescent - e.g. for cephalopods */
                temprec = bm->tot_SSB[species] * stock_prop[species][stock_id] * exp(recSTOCK[species][stock_id] * Ralpha_sp * (1.0 - bm->tot_SSB[species] * stock_prop[species][stock_id] / Rbeta_sp));

                /*
                if (do_debug && (bm->which_check == species)) {
                    fprintf( llogfp, "Time: %e, box%d-%d, species %s, ngene: %d, stock: %d, temprec: %e, totfishpop: %e, stock_prop: %e, recSTOCK: %e, Ralpha: %e, Larvae: %e, Rbeta: %e\n",
                            bm->dayt, bm->current_box, bm->current_layer, FunctGroupArray[species].groupCode,
                            ngene, stock_id, temprec, bm->tot_SSB[species], stock_prop[species][stock_id], recSTOCK[species][stock_id], Ralpha_sp,
                            EMBRYO[species].Larvae[stock_id][ngene][expect_id], Rbeta_sp);
                }
                 */
            break;
		default:
			quit("No such flagrecruit defined for vertebrates (%d) - value must be between 0 and 10 currently\n", recruit_sp);
			break;
		}

        // Check for any larval processes
        if ( FunctGroupArray[species].speciesParams[intersp_depend_recruit_id] > 0) {
            larval_scalar = Larval_Mortality(bm, species, stock_id, llogfp);
            temprec *= larval_scalar;
        }
        
        // Store the final number of recruits
		EMBRYO[species].BulkRecruits[ngene] = temprec;
        
        /*
        if (do_debug && (bm->which_check == species)) {
			fprintf( llogfp, "Time: %e, box%d-%d, species %s, ngene: %d, recruit_sp case: %d, BulkRecruits: %e, temprec: %e, Larvae: %e qid: %d)\n",
					bm->dayt, bm->current_box, bm->current_layer, FunctGroupArray[species].groupCode, ngene, recruit_sp,
					EMBRYO[species].BulkRecruits[ngene], temprec, EMBRYO[species].Larvae[stock_id][ngene][qid], qid);
		}
         */
	}
    
	return;
}
