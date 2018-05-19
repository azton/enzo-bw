// Parallel grid decomposition
 
// Robert Harkness
// April 2008
 
#include <mpi.h>
#include <hdf5.h>
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
 
#define MAXCPU 16384
#define FAIL 1
#define MAX_LINE_LENGTH 128
#define MAX_GRID_TAG_SIZE         16
#define MAX_TASK_TAG_SIZE         16
#define GRID_TAG_FORMAT        "8.8"
#define TASK_TAG_FORMAT        "8.8"
 
#include "macros_and_parameters.h"
 
// HDF5 prototypes
 
#include "extern_hdf5.h"
 
int Enzo_Dims_create(int nnodes, int ndims, int *dims);
 
 
Eint32 main(Eint32 argc, char *argv[])
{
 
  hid_t       file_id, dset_id, attr_id;
  hid_t       mem_dsp_id, file_dsp_id, attr_dsp_id;
  hid_t       file_type_id, mem_type_id;
  hid_t       int_file_type_id, int_mem_type_id;
  hid_t       file_acc_template;
  hid_t       xfer_prop_list;
 
  hsize_t     dims[4];
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  hsize_t     mem_stride, mem_count, attr_count;
  hsize_t     slab_stride[4], slab_count[4];
  hsize_t     slab_rank;
  hsize_t     slab_dims[4];
 
  hssize_t    mem_offset;
  hssize_t    slab_offset[4];
 
 
  int i, j, k, m, n;
  int a, b, c;
  int ii, jj, kk;
  int fln, ext, lex;
 
  int dim, rank, ngrids, gridcounter;
  int xpos1, ypos1, zpos1;
  int xpos2, ypos2, zpos2;
  int nx, ny, nz;

  int mode;
  char *argin;
  char s1[5];
  char proc_name[MAX_LINE_LENGTH];
 
  int enzo_layout[3] = {0,0,0};
  int out_dims[3];

  int Rank;
  int CompRank;
  int CompSize;
  int OutRank;
  int OutSize;
  int OutDims[3];
  int Dims[3]; 
  int Starts[3];
  int Ends[3];
  int TopGridDims[3];
 
  double DomainLeftEdge[3] = {0.0, 0.0, 0.0};
  double DomainRightEdge[3] = {1.0, 1.0, 1.0};
  double SubDomainLeftEdge[3];
  double SubDomainRightEdge[3];
  double Left[3];
  double Right[3];
  double x0, y0, z0, dx, dy, dz;
 
  double GridLeft[MAXCPU][3];
  double GridRight[MAXCPU][3];

  double *buff;
  double *obuff;
 
  char *GDin;
  char *GVin;
  char *Extension;
 
  char pid[MAX_GRID_TAG_SIZE];

  int io_log = 1;
  int io_log_d = 0;
 
//
 
  MPI_Arg mpi_argc;
  MPI_Arg mpi_size;
  MPI_Arg mpi_rank;
  MPI_Arg lname;

  int mpi_layout[3] = {0,0,0};
  int ncpu, jcpu;

  FILE *log;

//
 
  mpi_argc = argc;
 
  MPI_Init(&mpi_argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Get_processor_name(proc_name, &lname);

  fprintf(stderr, "Proc %d of %d is %s\n", mpi_rank, mpi_size, proc_name);

  ncpu = mpi_size;
  jcpu = mpi_rank;

/*
  Mode = 1 then Density and Velocity only
  Mode = 2 then Density only
  Mode = 3 then Velocity only

*/

  argin = argv[1];
  mode = 0;

  if (argc == 4) {
    strcpy(s1, "dv");
    i = strcmp(argin, s1);
    if ( i == 0 ) {
      mode = 1;
      GDin = argv[2];
      GVin = argv[3];
      if ( mpi_rank == 0 ) {
        printf("Input arg %s should be dv\n", argin);
        printf("GDin = %s\n", GDin);
        printf("GVin = %s\n", GVin);
        printf("Read Grid Density\n");
        printf("Read Grid Velocity\n");
      }
    }
  }

  if (argc == 4) {
    strcpy(s1, "d");
    i = strcmp(argin, s1);
    if ( i == 0 ) {
      mode = 2;
      GDin = argv[2];
      GVin = argv[3];
      if ( mpi_rank == 0 ) {
        printf("Input arg %s should be d\n", argin);
        printf("GDin = %s\n", GDin);
        printf("Read Grid Density\n");
      }
    }
  }

  if (argc == 4) {
    strcpy(s1, "v");
    i = strcmp(argin, s1);
    if ( i == 0 ) {
      mode = 3;
      GDin = argv[2];
      GVin = argv[3];
      if ( mpi_rank == 0 ) {
        printf("Input arg %s should be v\n", argin);
        printf("GVin = %s\n", GVin);
        printf("Read Grid Velocity\n");
      }
    }
  }

  if (mode == 0) {
    if ( mpi_rank == 0 ) {
      printf("Error: no matching keys\n");
    }
  }

 
  lex = 0;

  fln = strlen(GDin);
  ext = strcspn(GDin, ".");

  if ( fln-ext > 0 )
  {
    lex = strlen(strstr(GDin, "."));
    Extension = new char[MAX_LINE_LENGTH];
    strcpy(Extension, strstr(GDin, "."));
    lex = strlen(Extension);
  }

  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, jcpu);
 
  char *logname = new char[MAX_LINE_LENGTH];
  strcpy(logname, "GRIDlog");
  strcat(logname, pid);
  if (lex > 0)
    strcat(logname, Extension);
 
  log = fopen(logname, "a");
 
  char *GD = new char[MAX_LINE_LENGTH];
  strcpy(GD, "GD");
  strcat(GD, pid);
  if (lex > 0)
    strcat(GD, Extension);
 
  char *GV = new char[MAX_LINE_LENGTH];
  strcpy(GV, "GV");
  strcat(GV, pid);
  if (lex > 0)
    strcat(GV, Extension);
 
  mem_type_id = HDF5_R8;
  file_type_id = HDF5_FILE_R8;




  // Get attributes from GridDensity
 
  file_id = H5Fopen(GDin, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (io_log) fprintf(log, "H5Fopen with Name = %s\n", GDin);
    assert( file_id != h5_error );
 
  dset_id = H5Dopen(file_id, GDin);
    if (io_log) fprintf(log, "H5Dopen with Name = %s\n", GDin);
    assert( dset_id != h5_error );
 
  attr_id = H5Aopen_name(dset_id, "Component_Rank");
    if (io_log) fprintf(log, "H5Aopen with Name = Component_Rank\n");
    assert( attr_id != h5_error );
 
  h5_status = H5Aread(attr_id, HDF5_INT, &CompRank);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  attr_id = H5Aopen_name(dset_id, "Component_Size");
   if (io_log) fprintf(log, "H5Aopen with Name = Component_Size\n");
   assert( attr_id != h5_error );
 
  h5_status = H5Aread(attr_id, HDF5_INT, &CompSize);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  attr_id = H5Aopen_name(dset_id, "Rank");
    if (io_log) fprintf(log, "H5Aopen with Name = Rank\n");
    assert( attr_id != h5_error );

  h5_status = H5Aread(attr_id, HDF5_INT, &Rank);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  attr_count = 3;
 
  attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
    assert( attr_dsp_id != h5_error );

  attr_id = H5Aopen_name(dset_id, "Dimensions");
   if (io_log) fprintf(log, "H5Aopen with Name = Dimensions\n");
   assert( attr_id != h5_error );

  h5_status = H5Aread(attr_id, HDF5_INT, Dims);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  attr_id = H5Aopen_name(dset_id, "TopGridStart");
   if (io_log) fprintf(log, "H5Aopen with Name = TopGridStart\n");
   assert( attr_id != h5_error );
 
  h5_status = H5Aread(attr_id, HDF5_INT, Starts);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  attr_id = H5Aopen_name(dset_id, "TopGridEnd");
   if (io_log) fprintf(log, "H5Aopen with Name = TopGridEnd\n");
   assert( attr_id != h5_error );
 
  h5_status = H5Aread(attr_id, HDF5_INT, Ends);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  attr_id = H5Aopen_name(dset_id, "TopGridDims");
   if (io_log) fprintf(log, "H5Aopen with Name = TopGridDims\n");
   assert( attr_id != h5_error );
 
  h5_status = H5Aread(attr_id, HDF5_INT, TopGridDims);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(attr_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );


  ngrids = 1;
  rank = Rank;
 
  Enzo_Dims_create(ncpu, rank, mpi_layout);

  for (dim = 0; dim < rank; dim++)
  {
    enzo_layout[dim] = mpi_layout[rank-1-dim];
    ngrids *= enzo_layout[dim];
  }

  if ( jcpu == 0 ) {
    fprintf(stderr, "NumberOfGrids = %"ISYM"\n", ngrids);
    fprintf(stderr, "ENZO_layout %"ISYM" x %"ISYM" x %"ISYM"\n", enzo_layout[0], enzo_layout[1], enzo_layout[2]);
  }

  fprintf(log, "NumberOfGrids = %"ISYM"\n", ngrids);
  fprintf(log, "ENZO_layout %"ISYM" %"ISYM" %"ISYM"\n", enzo_layout[0], enzo_layout[1], enzo_layout[2]);

  fprintf(log, "Component_Rank = %"ISYM"\n", CompRank);
  fprintf(log, "Component_Size = %"ISYM"\n", CompSize);
  fprintf(log, "Rank = %"ISYM"\n", Rank); 
  fprintf(log, "Dims = %"ISYM" %"ISYM" %"ISYM"\n", Dims[0], Dims[1], Dims[2]);
  fprintf(log, "TopGridDims = %"ISYM" %"ISYM" %"ISYM"\n", TopGridDims[0], TopGridDims[1], TopGridDims[2]);
  fprintf(log, "TopGridStart = %"ISYM" %"ISYM" %"ISYM"\n", Starts[0], Starts[1], Starts[2]);
  fprintf(log, "TopGridEnd = %"ISYM" %"ISYM" %"ISYM"\n", Ends[0], Ends[1], Ends[2]);
 
  for (dim = 0; dim < rank; dim++)
  {
    SubDomainLeftEdge[dim] = Starts[dim] * (DomainRightEdge[dim]-DomainLeftEdge[dim])/((double) TopGridDims[dim]);
    SubDomainRightEdge[dim] = (Ends[dim]+1) * (DomainRightEdge[dim]-DomainLeftEdge[dim])/((double) TopGridDims[dim]);
    //  SubCellWidth[dim] = (SubDomainRightEdge[dim]-SubDomainLeftEdge[dim])/((double) enzo_layout[dim]);
  }
 
  x0 = SubDomainLeftEdge[0];
  y0 = SubDomainLeftEdge[1];
  z0 = SubDomainLeftEdge[2];
 
  dx = (SubDomainRightEdge[0]-SubDomainLeftEdge[0])/((double) enzo_layout[0]);
  dy = (SubDomainRightEdge[1]-SubDomainLeftEdge[1])/((double) enzo_layout[1]);
  dz = (SubDomainRightEdge[2]-SubDomainLeftEdge[2])/((double) enzo_layout[2]);
 
  a=enzo_layout[0];
  b=enzo_layout[1];
  c=enzo_layout[2];
 
  gridcounter = 0;
 
  for (kk = 0; kk < enzo_layout[2]; kk++)
    for (jj = 0; jj < enzo_layout[1]; jj++)
      for (ii = 0; ii < enzo_layout[0]; ii++)
      {
        n=gridcounter;
 
        if ( n == jcpu )
        {
      // rank to coordinate
        m = n;
        i = m/(b*c);
        m = m%(b*c);
        j = m/c;
        m = m%c;
        k = m;
      // coordinate to rank check
        m = ((i*b*c) + j*c) + k;
        fprintf(log,"Grid %"ISYM"  {%"ISYM" %"ISYM" %"ISYM"}  %"ISYM"\n",n,i,j,k,m);
 
        Left[0] =  x0 + dx * (double) ii;
        Right[0] = x0 + dx * (double) (ii+1);
        Left[1] =  y0 + dy * (double) jj;
        Right[1] = y0 + dy * (double) (jj+1);
        Left[2] =  z0 + dz * (double) kk;
        Right[2] = z0 + dz * (double) (kk+1);
 
        GridLeft[jcpu][0] = Left[0];
        GridLeft[jcpu][1] = Left[1];
        GridLeft[jcpu][2] = Left[2];
 
        GridRight[jcpu][0] = Right[0];
        GridRight[jcpu][1] = Right[1];
        GridRight[jcpu][2] = Right[2];
 
        for (m = 0; m < rank; m++)
          fprintf(log, "Grid %"ISYM"    Left   %16.8"FSYM"   Right  %16.8"FSYM"\n", gridcounter, Left[m], Right[m]);
 
        }
 
        gridcounter++;
 
      }


    xpos1 = nint(Left[0] * ((double) Dims[0]));
    ypos1 = nint(Left[1] * ((double) Dims[1]));
    zpos1 = nint(Left[2] * ((double) Dims[2]));

    xpos2 = nint(Right[0] * ((double) Dims[0])) - 1;
    ypos2 = nint(Right[1] * ((double) Dims[1])) - 1;
    zpos2 = nint(Right[2] * ((double) Dims[2])) - 1;

    nx = xpos2 - xpos1 + 1;
    ny = ypos2 - ypos1 + 1;
    nz = zpos2 - zpos1 + 1;

    fprintf(log, "Grid X: %"ISYM"  [%"ISYM",%"ISYM"]\n", nx, xpos1, xpos2);
    fprintf(log, "Grid Y: %"ISYM"  [%"ISYM",%"ISYM"]\n", ny, ypos1, ypos2);
    fprintf(log, "Grid Z: %"ISYM"  [%"ISYM",%"ISYM"]\n", nz, zpos1, zpos2);

    MPI_Barrier(MPI_COMM_WORLD);

    buff = new double[nx*ny*nz];

    obuff = new double[nx*ny*nz];
 
    dims[0] = CompRank;
    dims[1] = Dims[0];
    dims[2] = Dims[1];
    dims[3] = Dims[2];

    mem_type_id = HDF5_R8;

  // GridDensity

  // if mode = 1 or 2

  if ((mode == 1) || (mode == 2)) {
 
  dim = 0;

    //  file_acc_template = H5Pcreate (H5P_FILE_ACCESS);
    //    if (io_log) fprintf(log, "H5Pcreate file_access_template\n");
    //    assert( file_acc_template != h5_error );
 
    //  h5_status = H5Pset_fapl_mpio(file_acc_template, MPI_COMM_WORLD, MPI_INFO_NULL);
    //    if (io_log) fprintf(log, "H5Pset_fapl_mpio: %"ISYM"\n", h5_status);
    //    assert( h5_status != h5_error );
 
    file_acc_template = H5P_DEFAULT;
      if (io_log) fprintf(log, "Default file_access_template\n");
 
    file_id = H5Fopen(GDin, H5F_ACC_RDONLY, file_acc_template);
      if (io_log) fprintf(log, "H5Fopen with Name = %s\n", GDin);
      assert( file_id != h5_error );
 
    dset_id = H5Dopen(file_id, GDin);
      if (io_log) fprintf(log, "H5Dopen with Name = %s\n", GDin);
      assert( dset_id != h5_error );
 
    mem_offset = 0;
    mem_count = nx*ny*nz;
    mem_stride = 1;
 
    mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", mem_dsp_id);
      assert( mem_dsp_id != h5_error );
 
    h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
      if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    file_dsp_id = H5Screate_simple((Eint32) 4, dims, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", file_dsp_id);
      assert( file_dsp_id != h5_error );
 
    slab_offset[0] = dim;
    slab_stride[0] = 1;
    slab_count[0] = 1;
 
    slab_offset[1] = zpos1;
    slab_stride[1] = 1;
    slab_count[1] = nz;

    slab_offset[2] = ypos1;
    slab_stride[2] = 1;
    slab_count[2] = ny;

    slab_offset[3] = xpos1;
    slab_stride[3] = 1;
    slab_count[3] = nx;
 
    h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
      if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    //  xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
    //    if (io_log) fprintf(log, "H5Pcreate xfer_prop_list\n");
    //    assert( xfer_prop_list != h5_error );
 
    //  h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
    //    if (io_log) fprintf(log, "H5Pset_dxpl_mpio: %"ISYM"\n", h5_status);
    //    assert( h5_status != h5_error );
 
    xfer_prop_list = H5P_DEFAULT;
      if (io_log) fprintf(log, "Default xfer_prop_list\n");
 
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, &buff[0]);
      assert( h5_status != h5_error );
 
    //  h5_status = H5Pclose(xfer_prop_list);
    //    if (io_log) fprintf(log, "H5Pclose: %"ISYM"\n", h5_status);
    //    assert( h5_status != h5_error );
 
    h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(mem_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Fclose(file_id);
      if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    //  h5_status = H5Pclose(file_acc_template);
    //    if (io_log) fprintf(log, "H5Pclose: %"ISYM"\n", h5_status);
    //    assert( h5_status != h5_error );
 

  // Write subgrids

/* 
  int counter = 0;
  for ( k = 0; k < nz; k++) {
    for ( j = 0; j < ny; j++) {
      for ( i = 0; i < nx; i++) {
        obuff[counter] = buff[counter];
        counter++;
      }
    }
  }
*/

  // fprintf(stderr, "DATA INPUT OK\n");

  MPI_Barrier(MPI_COMM_WORLD);

  // Data in memory is considered 1D, stride 1, with zero offset
 
  mem_stride = 1;                   // contiguous elements
  mem_count = nx * ny * nz;         // number of elements in field
  mem_offset = 0;                   // zero offset in buffer
 
  // 1D memory model
 
  mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", mem_dsp_id);
    assert( mem_dsp_id != h5_error );
 
  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  OutRank = CompRank;
  OutSize = nx * ny * nz;
  OutDims[0] = nx;
  OutDims[1] = ny;
  OutDims[2] = nz;


  slab_rank = 4;
  slab_dims[0] = 1;
  slab_dims[1] = nx;
  slab_dims[2] = ny;
  slab_dims[3] = nz;

  slab_stride[0] = 1;      // contiguous elements
  slab_count[0] = 1;       // one component per call
  slab_offset[0] = 0;      // component Part of Npart
 
  slab_stride[1] = 1;      // contiguous elements
  slab_count[1] = nx;      // field dimensions
  slab_offset[1] = 0;      // complete field, no offset

  slab_stride[2] = 1;      // contiguous elements
  slab_count[2] = ny;      // field dimensions
  slab_offset[2] = 0;      // complete field, no offset

  slab_stride[3] = 1;      // contiguous elements
  slab_count[3] = nz;      // field dimensions
  slab_offset[3] = 0;      // complete field, no offset
 
  file_dsp_id = H5Screate_simple((Eint32) slab_rank, slab_dims, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", file_dsp_id);
    assert( file_dsp_id != h5_error );
 
  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

 
  if ( dim == 0 )
  {

    file_id = H5Fcreate(GD, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fcreate with name = %s\n", GD);
      assert( file_id != h5_error );
 
    dset_id =  H5Dcreate(file_id, GD, file_type_id, file_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Dcreate with Name = %s\n", GD);
      assert( dset_id != h5_error );
 
    attr_count = 1;
 
    attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "Component_Rank", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = Component_Rank\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_INT, &OutRank);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_id = H5Acreate(dset_id, "Component_Size", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = Component_Size\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_INT, &OutSize);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_id = H5Acreate(dset_id, "Rank", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = Rank\n");
      assert( attr_id != h5_error );

    h5_status = H5Awrite(attr_id, HDF5_INT, &Rank);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    attr_count = 3;
 
    attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );

    attr_id = H5Acreate(dset_id, "Dimensions", HDF5_FILE_I8, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = Dimensions\n");
      assert( attr_id != h5_error );

    h5_status = H5Awrite(attr_id, HDF5_INT, OutDims);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_id = H5Acreate(dset_id, "GridLeft", HDF5_FILE_R8, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = GridLeft\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_R8, Left);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_id = H5Acreate(dset_id, "GridRight", HDF5_FILE_R8, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = GridRight\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_R8, Right);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

  }
  else
  {
    file_id = H5Fopen(GD, H5F_ACC_RDWR, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fopen with Name = %s\n", GD);
      assert( file_id != h5_error );
 
    dset_id =  H5Dopen(file_id, GD);
      if (io_log) fprintf(log, "H5Dopen with Name = %s\n", GD);
      assert( dset_id != h5_error );
  }
 
  h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, buff);
    if (io_log) fprintf(log, "H5Dwrite: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(mem_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  } // mode = 1 or 2 

  //  GridVelocities

  //  if mode = 1 or 3

  if ((mode == 1) || (mode == 3)) {

  //  Get attributes from GridVelocities
 
  file_id = H5Fopen(GVin, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (io_log) fprintf(log, "H5Fopen with Name = %s\n", GVin);
    assert( file_id != h5_error );
 
  dset_id = H5Dopen(file_id, GVin);
    if (io_log) fprintf(log, "H5Dopen with Name = %s\n", GVin);
    assert( dset_id != h5_error );
 
  attr_id = H5Aopen_name(dset_id, "Component_Rank");
    if (io_log) fprintf(log, "H5Aopen with Name = Component_Rank\n");
    assert( attr_id != h5_error );
 
  h5_status = H5Aread(attr_id, HDF5_INT, &CompRank);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  attr_id = H5Aopen_name(dset_id, "Component_Size");
   if (io_log) fprintf(log, "H5Aopen with Name = Component_Size\n");
   assert( attr_id != h5_error );
 
  h5_status = H5Aread(attr_id, HDF5_INT, &CompSize);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  attr_id = H5Aopen_name(dset_id, "Rank");
    if (io_log) fprintf(log, "H5Aopen with Name = Rank\n");
    assert( attr_id != h5_error );

  h5_status = H5Aread(attr_id, HDF5_INT, &Rank);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  attr_count = 3;
 
  attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
    assert( attr_dsp_id != h5_error );

  attr_id = H5Aopen_name(dset_id, "Dimensions");
   if (io_log) fprintf(log, "H5Aopen with Name = Dimensions\n");
   assert( attr_id != h5_error );

  h5_status = H5Aread(attr_id, HDF5_INT, Dims);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  attr_id = H5Aopen_name(dset_id, "TopGridStart");
   if (io_log) fprintf(log, "H5Aopen with Name = TopGridStart\n");
   assert( attr_id != h5_error );
 
  h5_status = H5Aread(attr_id, HDF5_INT, Starts);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  attr_id = H5Aopen_name(dset_id, "TopGridEnd");
   if (io_log) fprintf(log, "H5Aopen with Name = TopGridEnd\n");
   assert( attr_id != h5_error );
 
  h5_status = H5Aread(attr_id, HDF5_INT, Ends);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  attr_id = H5Aopen_name(dset_id, "TopGridDims");
   if (io_log) fprintf(log, "H5Aopen with Name = TopGridDims\n");
   assert( attr_id != h5_error );
 
  h5_status = H5Aread(attr_id, HDF5_INT, TopGridDims);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(attr_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  fprintf(log, "Component_Rank = %"ISYM"\n", CompRank);
  fprintf(log, "Component_Size = %"ISYM"\n", CompSize);
  fprintf(log, "Rank = %"ISYM"\n", Rank);
  fprintf(log, "Dims = %"ISYM" %"ISYM" %"ISYM"\n", Dims[0], Dims[1], Dims[2]);
  fprintf(log, "TopGridDims = %"ISYM" %"ISYM" %"ISYM"\n", TopGridDims[0], TopGridDims[1], TopGridDims[2]);
  fprintf(log, "TopGridStart = %"ISYM" %"ISYM" %"ISYM"\n", Starts[0], Starts[1], Starts[2]);
  fprintf(log, "TopGridEnd = %"ISYM" %"ISYM" %"ISYM"\n", Ends[0], Ends[1], Ends[2]);


  dims[0] = CompRank;
  dims[1] = Dims[0];
  dims[2] = Dims[1];
  dims[3] = Dims[2];


  for ( dim = 0; dim < 3; dim++) {

    //  file_acc_template = H5Pcreate (H5P_FILE_ACCESS);
    //    if (io_log) fprintf(log, "H5Pcreate file_access_template\n");
    //    assert( file_acc_template != h5_error );

    //  h5_status = H5Pset_fapl_mpio(file_acc_template, MPI_COMM_WORLD, MPI_INFO_NULL);
    //    if (io_log) fprintf(log, "H5Pset_fapl_mpio: %"ISYM"\n", h5_status);
    //    assert( h5_status != h5_error );

    file_acc_template = H5P_DEFAULT;
      if (io_log) fprintf(log, "Default file_access_template\n");

    file_id = H5Fopen(GVin, H5F_ACC_RDONLY, file_acc_template);
      if (io_log) fprintf(log, "H5Fopen with Name = %s\n", GVin);
      assert( file_id != h5_error );
 
    dset_id = H5Dopen(file_id, GVin);
      if (io_log) fprintf(log, "H5Dopen with Name = %s\n", GVin);
      assert( dset_id != h5_error );
 
    mem_offset = 0;
    mem_count = nx*ny*nz;
    mem_stride = 1;
 
    mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", mem_dsp_id);
      assert( mem_dsp_id != h5_error );
 
    h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
      if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    file_dsp_id = H5Screate_simple((Eint32) 4, dims, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", file_dsp_id);
      assert( file_dsp_id != h5_error );

    slab_offset[0] = dim;
    slab_stride[0] = 1;
    slab_count[0] = 1;
 
    slab_offset[1] = zpos1;
    slab_stride[1] = 1;
    slab_count[1] = nz;

    slab_offset[2] = ypos1;
    slab_stride[2] = 1;
    slab_count[2] = ny;

    slab_offset[3] = xpos1;
    slab_stride[3] = 1;
    slab_count[3] = nx;
 
    h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
      if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    //  xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
    //    if (io_log) fprintf(log, "H5Pcreate xfer_prop_list\n");
    //    assert( xfer_prop_list != h5_error );
 
    //  h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
    //    if (io_log) fprintf(log, "H5Pset_dxpl_mpio: %"ISYM"\n", h5_status);
    //    assert( h5_status != h5_error );
 
    xfer_prop_list = H5P_DEFAULT;
      if (io_log) fprintf(log, "Default xfer_prop_list\n");
 
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, &buff[0]);
      assert( h5_status != h5_error );
 
    //  h5_status = H5Pclose(xfer_prop_list);
    //    if (io_log) fprintf(log, "H5Pclose: %"ISYM"\n", h5_status);
    //    assert( h5_status != h5_error );
 
    h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(mem_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Fclose(file_id);
      if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    //  h5_status = H5Pclose(file_acc_template);
    //    if (io_log) fprintf(log, "H5Pclose: %"ISYM"\n", h5_status);
    //    assert( h5_status != h5_error );
 

  // Write subgrids

/*
  int counter = 0;
  for ( k = 0; k < nz; k++) {
    for ( j = 0; j < ny; j++) {
      for ( i = 0; i < nx; i++) {
        obuff[counter] = buff[counter];
      }
    }
  }
*/

  // fprintf(stderr, "DATA INPUT2 OK\n");

  MPI_Barrier(MPI_COMM_WORLD);

  // Data in memory is considered 1D, stride 1, with zero offset
 
  mem_stride = 1;                   // contiguous elements
  mem_count = nx * ny * nz;         // number of elements in field
  mem_offset = 0;                   // zero offset in buffer
 
  // 1D memory model
 
  mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", mem_dsp_id);
    assert( mem_dsp_id != h5_error );
 
  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  OutRank = CompRank;
  OutSize = nx * ny * nz;
  OutDims[0] = nx;
  OutDims[1] = ny;
  OutDims[2] = nz;

  slab_rank = 4;
  slab_dims[0] = 3;
  slab_dims[1] = nx;
  slab_dims[2] = ny;
  slab_dims[3] = nz;

  slab_stride[0] = 1;      // contiguous elements
  slab_count[0] = 1;       // one component per call
  slab_offset[0] = dim;      // component Part of Npart
 
  slab_stride[1] = 1;      // contiguous elements
  slab_count[1] = nx;      // field dimensions
  slab_offset[1] = 0;      // complete field, no offset

  slab_stride[2] = 1;      // contiguous elements
  slab_count[2] = ny;      // field dimensions
  slab_offset[2] = 0;      // complete field, no offset

  slab_stride[3] = 1;      // contiguous elements
  slab_count[3] = nz;      // field dimensions
  slab_offset[3] = 0;      // complete field, no offset
 
  file_dsp_id = H5Screate_simple((Eint32) slab_rank, slab_dims, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", file_dsp_id);
    assert( file_dsp_id != h5_error );
 
  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  if ( dim == 0 )
  {

    file_id = H5Fcreate(GV, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fcreate with name = %s\n", GV);
      assert( file_id != h5_error );
 
    dset_id =  H5Dcreate(file_id, GV, file_type_id, file_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Dcreate with Name = %s\n", GV);
      assert( dset_id != h5_error );
 
    attr_count = 1;
 
    attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "Component_Rank", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = Component_Rank\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_INT, &OutRank);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_id = H5Acreate(dset_id, "Component_Size", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = Component_Size\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_INT, &OutSize);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    attr_id = H5Acreate(dset_id, "Rank", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = Rank\n");
      assert( attr_id != h5_error );

    h5_status = H5Awrite(attr_id, HDF5_INT, &Rank);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_count = 3;
 
    attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );

    attr_id = H5Acreate(dset_id, "Dimensions", HDF5_FILE_I8, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = Dimensions\n");
      assert( attr_id != h5_error );

    h5_status = H5Awrite(attr_id, HDF5_INT, OutDims);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_id = H5Acreate(dset_id, "GridLeft", HDF5_FILE_R8, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = GridLeft\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_R8, Left);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_id = H5Acreate(dset_id, "GridRight", HDF5_FILE_R8, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = GridRight\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_R8, Right);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

  }
  else
  {
    file_id = H5Fopen(GV, H5F_ACC_RDWR, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fopen with Name = %s\n", GV);
      assert( file_id != h5_error );
 
    dset_id =  H5Dopen(file_id, GV);
      if (io_log) fprintf(log, "H5Dopen with Name = %s\n", GV);
      assert( dset_id != h5_error );
  }
 
  h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, buff);
    if (io_log) fprintf(log, "H5Dwrite: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(mem_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  } // end of loop over dims

  } // if mode = 1 or 3
 
  fclose(log);
 
  MPI_Barrier(MPI_COMM_WORLD);

  if ( mpi_rank == 0 )
  {
    fprintf(stdout, "Grid decomposition complete\n");
  }

  MPI_Finalize();
 
}
