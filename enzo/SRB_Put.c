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

#define DATATYPE "garbage"
#define RESOURCE "sfs-tape-tgd"

#define BUFSIZE 4096*4096

#include "extern_srb.h"

/*
extern int clStatus(srbConn* conn);
extern void clFinish(srbConn* conn);

extern srbConn* srbConnect(char *srbHost, char* srbPort, char* srbAuth, char *userName,
                           char *domainName, char *authScheme,  char *serverDn);

extern int srbObjCreate(srbConn* conn, int catType, char *objID,
                        char *dataTypeName, char *resourceName, char *collectionName,
                        char *pathName, srb_long_t dataSize);

extern int srbObjClose(srbConn* conn, int desc);
extern int srbObjWrite(srbConn* conn, int desc, char* buf, int len);

extern void srb_perror(int fd, int error_id, char* error_mnemonic, int flags);
*/



int SRB_Put(char *srb_directory, char *srb_file, char *local_file)
{
    srbConn *conn;
    srbResult *res;
    int  bytesRead, bytesWritten, total=0, in_fd, out_fd;
    char *buf, cpFileName[512];
    long filebytes ;
    struct srbStat statbuf ;

    buf = malloc(BUFSIZE);

    conn = srbConnect (HOST, PORT, SRB_AUTH, USER, DOMAIN, AUTH_SCHEME, SERVER_DN);
    if (clStatus(conn) != CLI_CONNECTION_OK) {
        fprintf(stderr, "Connection to srbMaster failed.\n");
        fprintf(stderr, "%s", clErrorMessage(conn));
        srb_perror (2, clStatus(conn), "", SRB_RCMD_ACTION|SRB_LONG_MSG);
        clFinish (conn);
	free (buf);
        exit (1);
    }

    strcpy(cpFileName,  srb_file);
    out_fd = srbObjCreate (conn, MDAS_CATALOG, cpFileName, DATATYPE, RESOURCE ,srb_directory, NULL, -1);

    if (out_fd < 0) { /* error */
        fprintf(stderr, "Can't create srb file \"%s\", status = %d\n",cpFileName, out_fd);
        fprintf(stderr, "%s", clErrorMessage(conn));
        clFinish (conn);
	free (buf);
	exit (1);
    }

    in_fd = open (local_file, O_RDONLY, 0);
    if (in_fd < 0) { /* error */
        fprintf(stderr, "Can't open local file\"%s\"\n", local_file);
        clFinish (conn);
	free (buf);
        exit (1);
    }

    /* Read from the local file and write to the just created srb file */

    total = 0;

    while ((bytesRead = read(in_fd, buf, BUFSIZE)) > 0) {

        bytesWritten = srbObjWrite(conn, out_fd, buf, bytesRead);
        if (bytesWritten < bytesRead) {
           fprintf(stderr, "Error: Read %d bytes, Wrote %d bytes.\n ", bytesRead, bytesWritten);
            clFinish (conn);
	    free (buf);
            exit (1);
        }
	total += bytesWritten;
    }

    srbObjClose (conn, out_fd);
    close (in_fd);
    free (buf);
    clFinish(conn);
    return (0);
}
