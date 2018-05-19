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

extern int srbObjStat(srbConn* conn, int catType, char *path, struct stat *statbuf);

extern int srbObjUnlink (srbConn* conn, char *objID, char *collectionName);

extern void srb_perror(int fd, int error_id, char *error_mnenomic, int flags);
*/



int SRB_Delete(char *srb_directory, char *srb_file_or_dir)
{
    srbConn *conn;
    srbResult *res;
    int status;
    struct stat statbuf;
    char mypath[NAME_LEN];


    conn = srbConnect (HOST, PORT, SRB_AUTH, USER, DOMAIN, AUTH_SCHEME, SERVER_DN);
    if (clStatus(conn) != CLI_CONNECTION_OK) {
        fprintf(stderr,"Connection to srbMaster failed.\n");
        fprintf(stderr,"%s",clErrorMessage(conn));
        srb_perror (2, clStatus(conn), "", SRB_RCMD_ACTION|SRB_LONG_MSG);
        clFinish (conn);
        return (1);
    }

    sprintf (mypath, "%s/%s", srb_directory, srb_file_or_dir);
    status = srbObjStat(conn, MDAS_CATALOG, mypath, &statbuf);

    if ( status > 0 ) {
      fprintf(stderr, "Object status check failed severely!\n");
      clFinish(conn);
      return(1);
    }

    if (status < 0) { /* object doesn't exist? */
        fprintf(stderr, "Can't stat srb file %s, status = %d\n", mypath, status);
        fprintf(stderr,"%s",clErrorMessage(conn));
        clFinish (conn);
        return (0);
    }

    if ( status == 0 ) {
      if ((statbuf.st_mode & S_IFREG) != 0) {     /* A file */
          fprintf(stderr, "%s is a file with size = %lld\n", mypath, statbuf.st_size);
          fprintf(stderr, "Delete it! %s\n", mypath);
          status =  srbObjUnlink (conn, srb_file_or_dir, srb_directory);
          fprintf(stderr, "Delete status: %d\n", status);
          clFinish(conn);
          return (0);
      } else {
          fprintf (stderr, "%s is a directory - abort!\n", mypath);
          clFinish(conn);
          return (-1);
      }
    }

}
