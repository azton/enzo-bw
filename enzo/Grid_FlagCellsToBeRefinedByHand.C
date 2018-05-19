
//
// FindCellsToBeRefinedByHand.
// Reads arbitrary number of regions from a file.
// File changes as time progresses.
//
// For each level, RefineByHand expects two files 
//  (rather, one file and one set of files.)
// 1.) RefineByList.<level>
// 2.) RefineByHandFile.<level>.<file_number>
//
// <level> is the file for a given level.
// <file_number> is a sequential number for that level.
//
// RefineByList.level example:
// AfterTime 0.2 RefineByHandFile.0.0
// AfterTime 0.3 RefineByHandFile.0.1
//
// When t>0.2, this code will open RefineByHandFile.0.0, and flag 
// the cells described within. This will persist until the next file is found.
// (at and beyond t=0.3.)  
//
// The list/time file must be called RefineByList.[0-inf]
// The grid edge filename here is arbitrary, but should include
// level and sequence index for clarity.
//
// RefineByHandFile.0.0 example:
// LeftEdge = 0.0 0.0 0.0
// RightEdge = 0.1 0.1 0.1
// LeftEdge = 0.2 0.2 0.2
// RightEdge = 0.3 0.3 0.3
//
// Cell number not yet coded.
//

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"


int grid::FlagCellsToBeRefinedByHand(int level){

  //Open meta file FlagByHandList.level
  //Find latest set by time.  
  //Open latest file
  //Loop over file lines, reading left/right zone pairs.
  
  char line[MAX_LINE_LENGTH];
  char filename[MAX_LINE_LENGTH];
  sprintf(filename, "RefineByList.%d",level);
  int dbg = 2;
  if ( dbg > 0 )
    fprintf(stderr,"RefineByHand filename %s\n",filename);
  FILE *fptr = fopen(filename,"r");

  if( !fptr){
    //Quietly exit.
    return 0;
  }
  //fprintf(stderr,"found file\n");
  
  float timeTmp=-1.0;
  char *tmpRefineFile = new char[MAX_LINE_LENGTH];
  char *refineFile = new char[MAX_LINE_LENGTH];
  char *DefaultFileName = new char[MAX_LINE_LENGTH];
  sprintf(refineFile,"%s", "NoFileDefined");
  sprintf(DefaultFileName,"%s", "NoFileDefined");


  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    if( dbg > 1 )
      fprintf(stderr,"RefineByHand: line %s", line);
    if( sscanf(line, "AfterTime %"FSYM" %s", &timeTmp, tmpRefineFile) ){
      if( dbg > 0 )
	fprintf(stderr,"RefineByHand: read time %f is time %f\n", timeTmp, this->Time);
      if( this->Time > timeTmp )
	strcpy(refineFile, tmpRefineFile);
    }
  }  
  fclose(fptr);


  if( strcmp( refineFile, DefaultFileName ) == 0 ){
    //No file found, quietly return.
    return 0;
  }
  fprintf(stderr,"Flag By Hand: Opening File %s\n", refineFile);
  int NumberOfFlaggedCells = 0;

  FILE *fptr2 = fopen(refineFile,"r");
  if( fptr2 == NULL )
    return 0;

  int dim,index,i,j,k;                             //loop variables.
  int leftEdgeRead = FALSE, rightEdgeRead=FALSE;   //Booleans to keep track of what line has been read.
  float refineLeft[3], refineRight[3];             //tmp values for each region.
  int refineStart[3]={0,0,0}, refineEnd[3]={0,0,0};//tmp values for refinement.

  
  //Read left edge, read right edge, then flag cells.
  //Repeat until the end of the file.
  while (fgets(line, MAX_LINE_LENGTH, fptr2) != NULL) {


    if( sscanf(line, "LeftEdge = %"FSYM" %"FSYM" %"FSYM, 
	       refineLeft,refineLeft+1,refineLeft+2 ) ){
      leftEdgeRead = TRUE;
    }
    if( leftEdgeRead == TRUE )
      if( sscanf(line, "RightEdge = %"FSYM" %"FSYM" %"FSYM, 
		 refineRight,refineRight+1,refineRight+2 ) ){
	rightEdgeRead = TRUE;
      }
    if( rightEdgeRead == TRUE ){
      leftEdgeRead = FALSE; rightEdgeRead=FALSE;
      
      for( dim=0;dim<GridRank;dim++){
	refineStart[dim] = 
	  nint( (refineLeft[dim] - GridLeftEdge[dim] ) / CellWidth[dim][0] ) + GridStartIndex[dim];
	refineStart[dim] = max( GridStartIndex[dim], refineStart[dim]);
	refineEnd[dim] = 
	  nint( (refineRight[dim] - GridLeftEdge[dim] )/ CellWidth[dim][0]-1 ) + GridStartIndex[dim];
	refineEnd[dim] = min( GridEndIndex[dim], refineEnd[dim] );

      }
      

      for(k=refineStart[2];k<=refineEnd[2];k++)
	for(j=refineStart[1];j<=refineEnd[1];j++)
	  for(i=refineStart[0];i<=refineEnd[0];i++){
	    index = i + GridDimension[0]* (j + GridDimension[1]*k);
	    FlaggingField[index]++;
	    NumberOfFlaggedCells++;
	  }

    }
    
  }
  if( dbg > 0 )
    fprintf(stderr,"RefineByHand: Flagged %d\n",NumberOfFlaggedCells);
  fclose(fptr2);
  delete refineFile;
  delete tmpRefineFile;
  return NumberOfFlaggedCells;


}
