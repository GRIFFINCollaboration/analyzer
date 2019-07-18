#include <stdio.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#include "midas.h"
#include "web_server.h"


extern ANA_MODULE griffin_module;
extern ANA_MODULE descant_module;

char *analyzer_name = "Analyzer"; /* The analyzer name (client name)   */
INT analyzer_loop_period = 0; /* analyzer_loop call interval[ms](0:disable) */
INT odb_size = DEFAULT_ODB_SIZE; /* default ODB size */

ANA_MODULE *griffin_trigger_module[] = { &griffin_module, NULL };
ANA_MODULE *descant_trigger_module[] = { &griffin_module, NULL };

BANK_LIST griffin_ana_trigger_bank_list[] = { /* online banks */
  {"GRF3", TID_DWORD, 256, NULL}, {""} ,
};
BANK_LIST descant_ana_trigger_bank_list[] = { /* online banks */
  {"CAEN", TID_DWORD, 256, NULL}, {""} ,
};

ANALYZE_REQUEST analyze_request[] = {
	{"Trigger",                  /* equipment name */
		{1,                         /* event ID */
			TRIGGER_ALL,               /* trigger mask */
			//GET_SOME,                  /* get some events */ /* Removed ?? */
			GET_NONBLOCKING,           /* get some events */ /* later midas's */
			"SYSTEM",                  /* event buffer */
			TRUE,                      /* enabled */
			"", "",},
		NULL,                       /* analyzer routine */
		griffin_trigger_module,             /* module list */
		griffin_ana_trigger_bank_list,      /* bank list */
		1000,                       /* RWNT buffer size */
		TRUE,                       /* Use tests for this event */
	},
	{"VX1730",                  /* equipment name */
		{3,                         /* event ID */
			TRIGGER_ALL,               /* trigger mask */
			//GET_SOME,                  /* get some events */ /* Removed ?? */
			GET_NONBLOCKING,           /* get some events */ /* later midas's */
			"SYSTEM",                  /* event buffer */
			TRUE,                      /* enabled */
			"", "",},
		NULL,                       /* analyzer routine */
		descant_trigger_module,             /* module list */
		descant_ana_trigger_bank_list,      /* bank list */
		1000,                       /* RWNT buffer size */
		TRUE,                       /* Use tests for this event */
	}, {""},
};

/*--------------------------------------------------------------------------*/

volatile int shutdown_webserver;
static pthread_t web_thread;
INT analyzer_init(){
	// Create new thread which launches the Web server
	int a1=1;
	pthread_create(&web_thread, NULL,(void* (*)(void*))web_server_main, &a1);

	return SUCCESS;
}

INT analyzer_exit(){
	printf("waiting for server to shutdown ...\n");
	shutdown_webserver = 1;
	pthread_join(web_thread, NULL);
	return CM_SUCCESS;
}
INT ana_begin_of_run(INT run_number, char *error){ return CM_SUCCESS; }
INT ana_end_of_run(INT run_number, char *error){ return CM_SUCCESS; }
INT ana_pause_run(INT run_number, char *error){ return CM_SUCCESS; }
INT ana_resume_run(INT run_number, char *error){ return CM_SUCCESS; }
INT analyzer_loop(){ return CM_SUCCESS; }
