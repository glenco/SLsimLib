    /*
    double *common_scan_blockH5(H5std_string filename
			       , H5std_string dataset_name
			       , hsize_t n_rows
			       , hsize_t n_cols
			       , hsize_t offset_rows
			       , hsize_t *read_rows
			       ); 

    */

    bool getFilename (FILE *pFile, char *filename) {
    
      int MAXSIZE = 0xFFF;
      char proclnk[0xFFF];
      int fno;
      ssize_t r;
      if (pFile != NULL)
      {   
        fno = fileno(pFile);
        sprintf(proclnk, "/proc/self/fd/%d", fno);
        r = readlink(proclnk, filename, MAXSIZE);
        if (r < 0)
        {   
            printf("failed to readlink\n");
            return 0;
        }   
        filename[r] = '\0';
        printf ("FILENAME: %s\n", filename);
    
        return true;
      }   

      return false;
    } 

  
  // added by Marcos Pellejero & Antonio Dorta***********
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                             
   * Copyright by The HDF Group.                                               *    
   * Copyright by the Board of Trustees of the University of Illinois.         *                                                                               
   * All rights reserved.                                                      *                                                                                
   *                                                                           *                                                                                
   * This file is part of HDF5.  The full HDF5 copyright notice, including     *                                                                                
   * terms governing use, modification, and redistribution, is contained in    *                                                                                
   * the COPYING file, which can be found at the root of the source code       *                                                                                
   * distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.  *                                                                                
   * If you do not have access to either file, you may request a copy from     *                                                                                
   * help@hdfgroup.org.                                                        *                                                                                
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


  /*                                                                                                                                                            
   *   This example shows how to read data from a chunked dataset.                                                                                              
   *   We will read from the file created by extend.cpp                                                                                                         
   */
  /*Modified by Marcos and Antonio 4th July 2017*/

  double *common_scan_blockH5(H5std_string filename, H5std_string dataset_name, hsize_t n_rows, hsize_t n_cols, hsize_t offset_rows, hsize_t *read_rows)
  {
    const int RANK = 2;
    double *data_out = NULL;
    *read_rows = 0;
    // Try block to detect exceptions raised by any of the calls inside it                                                                                      
    try
      {
	/*                                                                                                                                                  
	 * Turn off the auto-printing when failure occurs so that we can                                                                                    
	 * handle the errors appropriately                                                                                                                  
	 */
	Exception::dontPrint();
	/*                                                                                                                                                  
	 * Open the file and the dataset.                                                                                                                   
	 */
	H5File file( filename, H5F_ACC_RDONLY );
	DataSet dataset = file.openDataSet(dataset_name);
	/*                                                                                                                                                  
	 * Get filespace for rank and dimension                                                                                                             
	 */
	DataSpace filespace = dataset.getSpace();
	/*                                                                                                                                                  
	 * Get number of dimensions in the file dataspace                                                                                                   
	 */
	int rank = filespace.getSimpleExtentNdims();
	/*                                                                                                                                                  
	 * Get and print the dimension sizes of the file dataspace                                                                                          
	 */
	hsize_t dims[RANK];    // dataset dimensions                                                                                                        
        hsize_t count_rows;
        hsize_t offset[RANK];
	rank = filespace.getSimpleExtentDims( dims );
	cout << "dataset rank = " << rank << ", dimensions "
	     << (hsize_t)(dims[0]) << " x "
	     << (hsize_t)(dims[1]) << endl;
	/*                                                                                                                                                  
	 * Define the memory space to read dataset.                                                                                                         
	 */
	hsize_t file_rows = dims[0], file_cols = dims[1];
	bool readALL = false;
	
	if (n_cols != file_cols) {
	  cout << "ERROR!!! You specified " << n_cols << " cols, but file contains " << dims[1] << " cols.\n";
	  return data_out;
	}
	
        if (offset_rows > file_rows) {
	  // Reading out of data                                                                                                                                  
	  cout << "ERROR!!! File has " << file_rows << " rows and you wanted to read till row " << offset_rows + n_rows << "\n";
	  return data_out;
	}
	else if ((offset_rows == 0) && (file_rows >= n_rows)) {
	  // Read ALL data                                                                                                                                        
	  count_rows = n_rows;
	  readALL = true;
	  cout << "Reading ALL\n";
	}
	else if (offset_rows + n_rows > file_rows) {
	  // We are out of bound, limit the max rows                                                                                                              
	  count_rows = file_rows - offset_rows;
	  cout << "Reading LIMIT ROWS: " << count_rows << " FROM ROW " << offset_rows << "\n";
	}
	else {
	  // Inside bounds                                                                                                                                        
	  count_rows = n_rows;
	  cout << "Reading INSIDE\n";
	}
	
	/*                                                                                                                                                  
	 * Read dataset back and display.                                                                                                                   
	 */
	
	
	// Allocate only the required memory
	data_out = new double [count_rows*file_cols]; 
	
	// READ DATA!!! (we are going to read count_rows x file_cols)
	if (readALL) {
	  DataSpace memspace(RANK, dims);
	  dataset.read(data_out, PredType::NATIVE_DOUBLE, memspace, filespace);
	}
	else {
	  hsize_t count[RANK];
	  count[0] = count_rows;
	  count[1] = file_cols;
	  offset[0] = offset_rows;
	  offset[1] = 0;
	  DataSpace memspace(RANK, count);
	  filespace.selectHyperslab(H5S_SELECT_SET, count, offset);
	  dataset.read(data_out, PredType::NATIVE_DOUBLE, memspace, filespace);
	}
	*read_rows = count_rows;
	
      }  // end of try block
    // catch failure caused by the H5File operations
    catch( FileIException error )
      {
	error.printError();
	*read_rows = 0;
	return NULL;
      }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
      {
	error.printError();
	*read_rows = 0;
	return NULL;
      }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
      {
	error.printError();
	*read_rows = 0;
	return NULL;
      }
    //return 0;
    return data_out;
  }
  

