#include <stdio.h>
#include <sys/file.h>
#include <sys/stat.h>
#include "srbClient.h"

#define HOST       	NULL
#define PORT            NULL
#define SRB_AUTH        NULL	/* For 'PASSWD_AUTH', 'ENCRYPT1' */
#define USER            NULL
#define DOMAIN          NULL
#define AUTH_SCHEME     NULL	/* 'PASSWD_AUTH', 'ENCRYPT1', 'GSI_AUTH' */
#define SERVER_DN       NULL

#define NAME_LEN 256

#include "extern_srb.h"

/*
extern int clStatus(srbConn* conn);
extern void clFinish(srbConn* conn);

extern srbConn* srbConnect(char *srbHost, char* srbPort, char* srbAuth, char *userName,
                           char *domainName, char *authScheme,  char *serverDn);

extern int srbCreateCollect(srbConn* conn, int catType, char* parentCollect, char* newCollect);

extern void srb_perror(int fd, int error_id, char* error_mnemonic, int flags);
*/



int SRB_Mkdir(char *srb_directory, char *srb_subdir)
{
    srbConn *conn;
    srbResult *res;
    struct stat statbuf;
    char mypath[NAME_LEN];
    int status;

    conn = srbConnect(HOST, PORT, SRB_AUTH, USER, DOMAIN, AUTH_SCHEME, SERVER_DN);
    if ( clStatus(conn) != CLI_CONNECTION_OK ) {
        fprintf(stderr, "Connection to srbMaster failed.\n");
        fprintf(stderr, "%s",clErrorMessage(conn));
        srb_perror(2, clStatus(conn), "", SRB_RCMD_ACTION|SRB_LONG_MSG);
        clFinish(conn);
        return (1);
    }

/* check if dir exists */

    sprintf (mypath, "%s/%s", srb_directory, srb_subdir);
    status = srbObjStat(conn, MDAS_CATALOG, mypath, &statbuf);

    if ( status > 0 ) {
      fprintf(stderr, "Object status check failed severely!\n");
      clFinish(conn);
      return(1);
    }

    if ( status == 0 ) { 
      if ((statbuf.st_mode & S_IFREG) != 0) {     /* A file */
          fprintf(stderr, "Target object %s is a file!\n", mypath);
          clFinish(conn);
          return (1);
      } else {
          fprintf(stderr, "Directory already exists!\n", mypath);
          clFinish(conn);
          return (0);
      }
    }

/* status < 0 == no object */

    if ( status < 0 ) {
      status = srbCreateCollect(conn, MDAS_CATALOG, srb_directory, srb_subdir);
      if ( status != 0 ) {
          fprintf(stderr, "SRB create directory %s failed\n", mypath);
          clFinish(conn);
          return (1);
      } else {
          fprintf(stderr, "New directory %s created\n", mypath);
          clFinish(conn);
          return (0);
      }
    }

}
