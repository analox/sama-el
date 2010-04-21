/*
    prepdoc1.cpp
    2002-05-07
    R. Alberts

    This file will prepare the doxygen documentation by changing the
    file "ArrayOp.h" of the Array package.
    This change is done by setting the defined operators
    "UnaryOperator", "BinaryOperator", "CompoundAssignmentOperator",
    "UnaryFunction" and "BinaryFunction" and their instances
    into comments. 
    This is necessary because doxygen is not able to deal with
    the kind of definition used for these operators and will
    create a corrupted reference for "ArrayOp.h".

 */


#include <string>
#include <fstream>
#include <iostream>

using namespace std;


fstream input_file,  // The current input file used.
        output_file; // The current output file used. 

const unsigned buf_size( 16384 );
char           buf[ buf_size ]; 


// Will copy the whole content of the source file with
// name "input_name" to the destination file with
// name "copy_name":
//
void make_copy( string input_name, string copy_name )
{
    ifstream       inputfile;
    ofstream       copyfile;

    inputfile.open( input_name.c_str( ) );
    if ( !inputfile )
    {
        cout << "Can't find file \"" << input_name << "\"!" << endl;
        exit( -1 );
    }

    inputfile.seekg( 0, ios::beg );    

    copyfile.open( copy_name.c_str( ), ios::trunc );
    if ( !copyfile )
    {
        cout << "Can't create file \"" << copy_name << "\"!" << endl;
        exit( -1 );
    }

    do {
        if ( inputfile.getline( buf, buf_size ) )
        {
	  copyfile << buf << endl;
        }
    } while ( inputfile );

    inputfile.close( );
    copyfile.close( );
}


int main()
{
    string   directory,
             filename,             
             copyfilename,
             fullname,
             starttoken,
             endtoken,
             buffer;               
    bool     out( false);


    directory    = "../../include/Array/";
    filename     = "ArrayOp.h";
    copyfilename = "./arraytmpfile";
    starttoken   = "//>> Begin of defined operators";
    endtoken     = "//>> End of defined operators";
    
    fullname = directory + filename;    
    make_copy( fullname, copyfilename );

    input_file.open( copyfilename.c_str( ), fstream::in );
    if ( !input_file )
    {
	cout << "Can't find file \"" << copyfilename << "\"!" << endl;
        exit( -1 );
    }
    output_file.open( fullname.c_str( ), fstream::out );
    if ( !output_file )
    {
	cout << "Can't find file \"" << fullname << "\"!" << endl;
        exit( -1 );
    }

    // The section with the defined operators is marked by
    // the two comments "//>> Begin of defined operators"
    // and "//>> End of defined operators" enclosing this
    // section. All lines of the file between these
    // to comments will be changed by placing the comment
    // marker "//" before them:
    do {
        if ( input_file.getline( buf, buf_size ) )
        {
            buffer = ( string ) buf;
            if ( buffer.find( starttoken, 0 ) != string::npos )
	    {
	        out = true;
	    } 
            else if ( buffer.find( endtoken, 0 ) != string::npos )
	    {
	        out = false;
	    }
            else 
	    {
	        if ( out )
	        {
		    buffer = "//" + buffer;                  
	        } 
                output_file << buffer << endl;
	    }
        }
    } while( !input_file.eof( ) );


    input_file.close( );
    output_file.close( );

    return 0;
}









