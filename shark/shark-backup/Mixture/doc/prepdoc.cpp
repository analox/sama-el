/*
    prepdoc.cpp
    2002-04-30
    R. Alberts

    This file will prepare the doxygen documentation by extracting the names
    of the doxygen files with the class/namespace/file descriptions and 
    creating the html-start-file "[package-name]_Reference.html" from file 
    "preindex.ht_" where placeholders are used to mark the positions where the
    extracted filenames has to be inserted.
    This is necessary because doxygen not always creates files with the
    same name.

 */



#include <string>
#include <fstream>
#include <iostream>

using namespace std;


// Name of the package:
const string  package_name   = "Mixture";


// Name of the reference index file:
const string  final_filename = "./" + package_name + "_Reference.html";


fstream curr_input_file,  // The current input file used.
        final_file;       // The final reference index file. 


// Total number of placeholders used:
const unsigned no_placeholders = 7;


// List of all placeholders - one placeholder in each line:
//--------------------------------------------------------------------
// 1. position in line: The name of the class/file/namespace
//                      that is connected with the placeholder
// 2. position in line: The placeholder name
// 3. position in line: The doxygen file, where it will be
//                      searched for the name given at the 1. position.
//                      If this filename begins with a "$",
//                      then the filename will be the one
//                      at the 4. position of the line
//                      with those name at the 1. position
//                      that follows the "$".
// 4. position in line: Initially empty. Here the final
//                      filename, that will replace the
//                      placeholder will be stored.
// 
string placeholder_list[ no_placeholders * 4 ] = {
  "CodeBook.cpp", "01", "files.html", "",
  "KernelDensity.cpp", "02", "files.html", "",
  "LocalLinearRegression.cpp", "03", "files.html", "",
  "MixtureOfGaussians.cpp", "04", "files.html", "",
  "RBFN-PTO.cpp", "05", "files.html", "",
  "fastem.cpp", "06", "files.html", "",
  "rbfn-example.cpp", "08", "examples.html", ""
};



// Searches in doxygen file that was opened before with
// the file handle "searchfile" for a
// string "searchname" and the corresponding link
// to the doxygen file which contains the description of
// the class. The name of this file is stored in "filename"
//
void searchFile( fstream &searchfile, string searchname, string &filename )
{
    bool              found( false );   // searchname found?
    string::size_type curr_pos;         // current position in buffer with
                                        // current line of searchfile
    char              buf[ 16384 ];     // buffer with current line of
                                        // searchfile ( C string )
    string            buffer;           // buffer with current line of
                                        // searchfile ( C++ string )
    char              start,            // character before the found 
                                        // searchname
                      end,              // character after the found 
                                        // searchname
                      curr_char;        // current character read
    unsigned          search_start_pos; // current position, where 
                                        // the search for the name in the
                                        // current line starts 


    // go to begin of searchfile:
    searchfile.seekg( 0, ios::beg );    

    do {
        // get the next line of the searchfile:
        if ( searchfile.getline( buf, 16384 ) )
        {
            buffer = ( string ) buf;

            // start search at beginning of the current line
            search_start_pos = 0;

            do {
                curr_pos = buffer.find( searchname, search_start_pos );

                if ( curr_pos != string::npos )
	        {
		    // If searchname occurrs more than one time:
                    search_start_pos = curr_pos + 1;

                    // Get characters before and after the found 
                    // searchname:
                    start = buffer[ curr_pos - 1 ];
                    end = buffer[ curr_pos + searchname.size( ) ];
                }

                // Found searchname must be enclosed by HTML link
		// command ( <a href="...">searchname['&' for
                // templates]</a> ):
	        if ( curr_pos != string::npos &&
                     start == '>' && ( end == '<' || end == '&' ) )
	        {
	            // seek start position of link to the file
	            // with the description of the class/namespace/file
	            curr_pos--;
	            do {
		        curr_pos--;
                        curr_char = buffer[ curr_pos ];
	            } while( curr_char != '=' );
		    // Go to first character of filename
                    do
		    {
		        curr_pos++;
                        curr_char = buffer[ curr_pos ];
		    } while ( curr_char != '\"' ); 
                    // copy characters of filename
                    do
		    {
		        curr_pos++;
                        curr_char = buffer[ curr_pos ];
                        if ( curr_char != '\"' )
		        {
                            filename += curr_char;
                        }
		    } while ( curr_char != '\"' ); 
		    found = true;
	        }
            } while ( curr_pos != string::npos );
        }
    } while ( !searchfile.eof() && !found );

    if ( !found )
    {
        cout << "Can't find link for class/namespace/file \"" << searchname 
             << "\"!" << endl;
	exit( -1 );
    }
}


// Given a buffer string containing the current line of a file,
// this function searches for the occurrences of the string
// "placeholder" and replaces it by "linkname":
//
void replace_link( string &line, string placeholder, string linkname )
{
    // current position of found "placeholder" in "line"
    string::size_type found_pos( 0 );

    // Replace all occurrences of "placeholder" in the "line" string:
    do {
        // Get next occurrence position:
        found_pos = line.find( placeholder, found_pos );
 
        // Found another occurrence? Then replace it:
        if ( found_pos != string::npos )
	{
	    line.replace( found_pos, placeholder.size( ), linkname ); 
	}
    } while ( found_pos != string::npos );

    return;
}


int main()
{

    string   directory,            // the directory, where the
                                   // doxygen files can be found
             name,                 // the name to search for
                                   // (1. position in the placeholder list) 
             placeholder,          // the name of the placeholder
                                   // (2. position in the placeholder list) 
             filename,             // the name of the file to search in
                                   // (3. position in the placeholder list) 
             finalname,            // the name of the file with the
                                   // description for the class/namespace/file
                                   // that will replace the placeholder
                                   // (4. position in the placeholder list) 
             nestedname;           // the name that is given as filename
                                   // with an introducing "$" (see placeholder
                                   // list)
    unsigned curr_placeholder_no,  // current (line) number of placeholder
                                   // list
             curr_phn2;            // second counter for
                                   // placeholders ((p)lace(h)older (n)o (2))
    char     buf[ 1024 ];          // buffer for current line a file 
                                   // ( C-style )
    string   buffer;               // buffer for current line a file 
                                   // ( C++-style )


    directory = "./html/";

    for ( curr_placeholder_no = 0; 
          curr_placeholder_no < no_placeholders;
          curr_placeholder_no++                  )
    {
        // Get name of the file to search in:
        filename = placeholder_list[ curr_placeholder_no * 4 + 2 ];
        if ( filename[ 0 ] == '$' )
	{
	    nestedname.assign( filename, 1, filename.size( ) - 1 );
	    for ( curr_phn2 = 0; curr_phn2 < curr_placeholder_no; curr_phn2++ )
	    {
	        if ( placeholder_list[ curr_phn2 * 4 ] == nestedname )
		{
		    filename = placeholder_list[ curr_phn2 * 4 + 3 ];
		}
            }
        }

        // Open file for search:
        filename = directory + filename;
        curr_input_file.open( filename.c_str( ), fstream::in );

        if ( !curr_input_file )
        {
	    cout << "Can't find file \"" << filename 
                 << "\" - maybe doxygen was not started yet?" << endl;
            exit( -1 );
        }

        // search for name:
        name = placeholder_list[ curr_placeholder_no * 4 ];

        searchFile( curr_input_file, name, 
                    placeholder_list[ curr_placeholder_no * 4 + 3 ] );

        curr_input_file.close( );
    }

    // open file where the extracted filenames shall be inserted
    // and the file that will contain the final version:
    curr_input_file.open( "./preindex.ht_", fstream::in );
    if ( !curr_input_file )
    {
        cout << "Can't find file \"./preindex.ht_\"!" << endl;
        exit( -1 );
    }    
    final_file.open( final_filename.c_str( ), fstream::out );
    if ( !final_file )
    {
        cout << "Can't create file \"" << final_filename << "\"!" << endl;
        exit( -1 );
    }    

    do {
        // get the next line of the file:
        if ( curr_input_file.getline( buf, 1024 ) )
        {
            buffer = ( string ) buf;

            // Replace all placeholders by the extracted file names:
            for ( curr_placeholder_no = 0; 
                  curr_placeholder_no < no_placeholders;
                  curr_placeholder_no++                  )
            {
                placeholder = "#" + 
                              placeholder_list[ curr_placeholder_no * 4 + 1 ];
                finalname = placeholder_list[ curr_placeholder_no * 4 + 3 ];

                replace_link( buffer, placeholder, finalname );

            }

            // Write line with replacements to the final file:
	    final_file << buffer << endl;
        }
    } while ( !curr_input_file.eof( ) );

    curr_input_file.close( );
    final_file.close( );

    return 0;
}


