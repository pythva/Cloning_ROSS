#include "clone_test_2.h"


static MPI_Datatype CLONE_STATE_TYPE = MPI_DATATYPE_NULL;
static MPI_Request  clone_rreq       = MPI_REQUEST_NULL;
//static int          clone_rreq_posted = 0;

void build_clone_state_type(void)
{
    /* struct is one int → contiguous */
    MPI_Type_contiguous(1, MPI_INT, &CLONE_STATE_TYPE);
    MPI_Type_commit(&CLONE_STATE_TYPE);
    //printf("hello worldddddddddddddddddddddddddddddddddddddd\n");
}


void clone_send_state(int receiver)
{
    //if (g_tw_mynode != 0) return;

    tw_lp *lp0 = g_tw_lp[0];
    /*
    if (!lp0)
        fprintf(stderr, "[rank %d] lp0 is NULL!\n", g_tw_mynode);
    else if (!lp0->cur_state)
        fprintf(stderr, "[rank %d] lp0->cur_state is NULL!\n", g_tw_mynode);
    else {
        clone_test_2_state *state =
            (clone_test_2_state *) lp0->cur_state;
        fprintf(stderr,
                "[rank %d] counter before send = %d (addr %p)\n",
                g_tw_mynode, state->counter, (void*)state);
    }*/
    clone_test_2_state *state =
        (clone_test_2_state *) lp0->cur_state;

    //MPI_Request sreq;
    //MPI_Isend(state, 1, CLONE_STATE_TYPE,
      //        /*dest=*/1, /*tag=*/200,
      //        MPI_COMM_WORLD, &sreq);
    //MPI_Wait(&sreq, MPI_STATUS_IGNORE);
    MPI_Send(state, 1, CLONE_STATE_TYPE,
             /*dest=*/receiver, /*tag=*/200, MPI_COMM_WORLD);
}

void clone_progress_recv(int sender)
{
    //if (g_tw_mynode != 1) return;      /* only rank 1 does the work */

    tw_lp *lp1 = g_tw_lp[0];
    /*
    if (!lp1)
        fprintf(stderr, "[rank %d] lp1 is NULL!\n", g_tw_mynode);
    else if (!lp1->cur_state)
        fprintf(stderr, "[rank %d] lp1->cur_state is NULL!\n", g_tw_mynode);
    else {
        clone_test_2_state *state =
            (clone_test_2_state *) lp1->cur_state;
        fprintf(stderr,
                "[rank %d] counter *before* recv = %d (addr %p)\n",
                g_tw_mynode, state->counter, (void*)state);
    }*/
    clone_test_2_state *state =
        (clone_test_2_state *) lp1->cur_state;

    /* Blocking until the message from rank 0 arrives */
    MPI_Recv(state, 1, CLONE_STATE_TYPE,
             /*src=*/sender, /*tag=*/200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

tw_peid
clone_test_2_map(tw_lpid gid)
{
  return (tw_peid) gid / g_tw_nlp;
}

void
clone_test_2_init(clone_test_2_state * s, tw_lp * lp)
{
  s->counter++; //all the lp will increment one 

  // now only the lp 0 will send a message to lp 0
  
  if (lp->gid == 0)
  {
   /*
	 tw_event_send(
    tw_event_new(lp->gid,
           tw_rand_exponential(lp->rng, DUNE_PILOT_START_TIME),
           lp));
  */
	  tw_event *e;
	  clone_test_2_message *m;
	  e = tw_event_new(lp->gid,
           tw_rand_exponential(lp->rng, DUNE_PILOT_START_TIME),
           lp);
	  m = tw_event_data(e);
	  m->type_of_message = regular;
	  tw_event_send(e);

	}
}

void
clone_test_2_event_handler(clone_test_2_state * s, tw_bf * bf, clone_test_2_message * m, tw_lp * lp)
{
  // assuming only lp 0 send message to lp0, then lp0 and lp0 is the only person that do the incrementation

  s->counter++; //only lp0 should the lp will increment one 

}

void
clone_clone_test_2_event_handler(clone_test_2_state * s, tw_bf * bf, clone_test_2_message * m, tw_lp * lp, int * clone_handler)
{ 
  // assuming only lp 0 send message to lp0, then lp0 and lp0 is the only person that do the incrementation
  /*
  if((lp->id) == 3){
  	printf("INSIDE LP3 HANDLER");
  }

  printf("PE ID: %d\n", g_tw_mynode);
  */

  /*
  if(((s->counter) % 2 == 1)){
  	*clone_handler = 1;  
  }
  s->counter++; //only lp0 should the lp will increment one 
  if((s->counter) < 3){//three branch will be 7, two branch is 5, if one branch then 3
  	tw_event_send(
    		tw_event_new(lp->gid,
           	tw_rand_exponential(lp->rng, DUNE_PILOT_START_TIME),
           	lp));
  }
*/	

	//tw_event *event;
        //clone_test_2_message *message;
 	
	switch(m->type_of_message)
    {
    	case branch:
    	{
		// current just incremnet
		// need to distinguish which is main with is differnet
		// now just send regular
		//s->counter++;// for now

		if(m->sender_pe == g_tw_mynode){// origional
			s->counter+=100;
		}else{// branc
			s->counter++;
		}
		if((s->counter) < 3){
		tw_event *event;
                clone_test_2_message *message;
                event = tw_event_new(lp->gid,
                tw_rand_exponential(lp->rng, DUNE_PILOT_START_TIME),
                                lp);
                message = tw_event_data(event);
                message->type_of_message = regular;
		message->sender_pe = g_tw_mynode;// current pe
                tw_event_send(event);}


	break;
	}
	case regular:
	{
		if(((s->counter) % 2 == 1)){// case where need to branch
        		*clone_handler = 1;
			// send with branh type message
			tw_event *event;
          		clone_test_2_message *message;
          		event = tw_event_new(lp->gid,
           		tw_rand_exponential(lp->rng, DUNE_PILOT_START_TIME),
           		lp);
          		message = tw_event_data(event);
          		message->type_of_message = branch;// branch type message 
          		message->sender_pe = g_tw_mynode;// current pe
			tw_event_send(event);
			break;
  		}

		s->counter++; //regular situation, just increment

		// regular send
		if((s->counter) < 3){
			tw_event *event;
                	clone_test_2_message *message;
                	event = tw_event_new(lp->gid,
                	tw_rand_exponential(lp->rng, DUNE_PILOT_START_TIME),
                        	lp);
                	message = tw_event_data(event);
                	message->type_of_message = regular;
			message->sender_pe = g_tw_mynode;// current pe
               		tw_event_send(event);
		}
	break;
	}	
    }
}

void
clone_test_2_event_handler_rc(clone_test_2_state * s, tw_bf * bf, clone_test_2_message * m, tw_lp * lp)
{
  /*************************************************************************/
  /* Self-Initiated Model - so it will never rollback for now **************/
  /*************************************************************************/
  tw_error(TW_LOC, "Trying to Rollback a Self-Initiated DUNE Model - need to implement RB now\n");
  return;
}

void clone_test_2_commit(clone_test_2_state * s, tw_bf * bf, clone_test_2_message * m, tw_lp * lp)
{
    (void) s;
    (void) bf;
    (void) m;
    (void) lp;
}

void
clone_test_2_finish(clone_test_2_state * s, tw_lp * lp)
{
 printf("at lp: %lu, the counter is %d, CURRENT TIME IS: %f\n",lp->gid,s->counter, tw_now(lp)); 
}

tw_lptype       mylps[] = {
  {(init_f) clone_test_2_init,
   (pre_run_f) NULL,
   (event_f) clone_test_2_event_handler,
   (clone_event_f) clone_clone_test_2_event_handler,
   (revent_f) clone_test_2_event_handler_rc,
   (commit_f) clone_test_2_commit,
   (final_f) clone_test_2_finish,
   (map_f) clone_test_2_map,
  sizeof(clone_test_2_state)},
  {0},
};

void event_trace(clone_test_2_message *m, tw_lp *lp, char *buffer, int *collect_flag)
{
    (void) m;
    (void) lp;
    (void) buffer;
    (void) collect_flag;
    return;
}

void clone_test_2_stats_collect(clone_test_2_state *s, tw_lp *lp, char *buffer)
{
    (void) s;
    (void) lp;
    (void) buffer;
    return;
}

st_model_types model_types[] = {
    {(ev_trace_f) event_trace,
     0,
    (model_stat_f) clone_test_2_stats_collect,
    sizeof(int),
    NULL, //(sample_event_f)
    NULL, //(sample_revent_f)
    0},
    {0}
};

const tw_optdef app_opt[] =
{
  TWOPT_GROUP("DUNE Model: Defaults to 39,500 jobs so arrange LPs accordingly"),
  TWOPT_DOUBLE("remote", percent_remote, "desired remote event rate"),
  TWOPT_UINT("nlp", nlp_per_pe, "number of LPs per processor"),
  TWOPT_DOUBLE("mean", mean, "exponential distribution mean for timestamps"),
  TWOPT_DOUBLE("mult", mult, "multiplier for event memory allocation"),
  TWOPT_DOUBLE("lookahead", lookahead, "lookahead for events"),
  TWOPT_UINT("start-events", g_clone_test_2_start_events, "number of initial messages per LP"),
  TWOPT_UINT("stagger", stagger, "Set to 1 to stagger event uniformly across 0 to end time."),
  TWOPT_UINT("memory", optimistic_memory, "additional memory buffers"),
  TWOPT_CHAR("run", run_id, "user supplied run name"),
  TWOPT_END()
};


int
main(int argc, char **argv)
{

#ifdef TEST_COMM_ROSS
    // Init outside of ROSS
    //MPI_Init(&argc, &argv);
    //build_clone_state_type();       /* must precede any send/recv */
    // Split COMM_WORLD in half even/odd
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm split_comm;
    MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%2, mpi_rank, &split_comm);
    if(mpi_rank%2 == 1){
        // tests should catch any MPI_COMM_WORLD collectives
        MPI_Finalize();
    }
    // Allows ROSS to function as normal
    tw_comm_set(split_comm);
#endif

  unsigned int i;
  //build_clone_state_type();       /* must precede any send/recv */
  // set a min lookahead of 1.0
  lookahead = 1.0;
  tw_opt_add(app_opt);
  tw_init(&argc, &argv);
  
#ifdef USE_DAMARIS
    if(g_st_ross_rank)
    { // only ross ranks should run code between here and tw_run()
#endif
  if( lookahead > 1.0 )
    tw_error(TW_LOC, "Lookahead > 1.0 .. needs to be less\n");

  //reset mean based on lookahead
        mean = mean - lookahead;

  offset_lpid = g_tw_mynode * nlp_per_pe;
  ttl_lps = tw_nnodes() * nlp_per_pe;
  g_tw_events_per_pe = (mult * nlp_per_pe * g_clone_test_2_start_events) +
        optimistic_memory;
  //g_tw_rng_default = TW_FALSE;
  g_tw_lookahead = lookahead;

  tw_define_lps(nlp_per_pe, sizeof(clone_test_2_message));

  for(i = 0; i < g_tw_nlp; i++)
    {
    tw_lp_settype(i, &mylps[0]);
        st_model_settype(i, &model_types[0]);
    }

        if( g_tw_mynode == 0 )
    {
      printf("========================================\n");
      printf("clone_test_2 Model Configuration..............\n");
      printf("   Lookahead..............%lf\n", lookahead);
      printf("   Start-events...........%u\n", g_clone_test_2_start_events);
      printf("   stagger................%u\n", stagger);
      printf("   Mean...................%lf\n", mean);
      printf("   Mult...................%lf\n", mult);
      printf("   Memory.................%u\n", optimistic_memory);
      printf("   Remote.................%lf\n", percent_remote);
      printf("========================================\n\n");
    }

  tw_run();
#ifdef USE_DAMARIS
    } // end if(g_st_ross_rank)
#endif
  //printf("END!!!!!!");
  tw_end();
  
  return 0;
}

