extern int clStatus(srbConn* conn);

extern void clFinish(srbConn* conn);

extern srbConn* srbConnect(char *srbHost, char* srbPort, char* srbAuth, char *userName,
                           char *domainName, char *authScheme,  char *serverDn);

extern int srbCreateCollect(srbConn* conn, int catType, char* parentCollect, char* newCollect);

extern int srbObjCreate(srbConn* conn, int catType, char *objID,
                        char *dataTypeName, char *resourceName, char *collectionName,
                        char *pathName, srb_long_t dataSize);

extern int srbObjWrite(srbConn* conn, int desc, char* buf, int len);

extern int srbObjClose(srbConn* conn, int desc);

extern int srbObjStat(srbConn* conn, int catType, char *path, struct stat *statbuf);

extern int srbObjUnlink (srbConn* conn, char *objID, char *collectionName);

extern void srb_perror(int fd, int error_id, char* error_mnemonic, int flags);



