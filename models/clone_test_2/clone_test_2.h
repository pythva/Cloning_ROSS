#ifndef INC_clone_test_2_h
#define INC_clone_test_2_h

#include <ross.h>

typedef struct clone_test_2_state clone_test_2_state;
typedef struct clone_test_2_message clone_test_2_message;
typedef enum message_type message_type;

struct clone_test_2_state
{

  int counter;          // simple counter
};

enum message_type
  {
    hello_0 = 1,
    hello_1, 
  };


struct clone_test_2_message
{
  message_type     type_of_message;
	
};

#define DUNE_PILOT_START_TIME 1.0

tw_stime lookahead = 1.0;
static unsigned int stagger = 0;
static unsigned int offset_lpid = 0;
static tw_stime mult = 1.4;
static tw_stime percent_remote = 0.25;
static unsigned int ttl_lps = 0;
static unsigned int nlp_per_pe = 8;
static unsigned int g_clone_test_2_start_events = 1;
static unsigned int optimistic_memory = 1024;

// rate for timestamp exponential distribution
static tw_stime mean = 1.0;

static char run_id[1024] = "DUNE/NOVA model";


#endif

