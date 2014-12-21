#include "vars.h"
#include<iostream>
#include<fstream>
#include<sstream>

inline std::string Trim(const std::string& str, const std::string& delimiters = " \f\n\r\t\v" )
{
    const size_t f = str.find_first_not_of( delimiters );
    return f == std::string::npos ?
                "" :
                str.substr( f, str.find_last_not_of( delimiters ) + 1 );
}

void ParseVarsFile(const string& filename, double* dx, double* dy,
	                int* nRows, int* nCols, string* model, bool* input_mode, string& input_name)
{
	ifstream f(filename.c_str());
    
    if( f.is_open() )
    {
        while( !f.bad() && !f.eof())
        {
            const int c = f.peek();
            
            if( isspace(c) )
            {
                // ignore leading whitespace
                f.get();
            }else{
                if( c == '#' || c == '%' )
                {
                    // ignore lines starting # or %
                    string comment;
                    getline(f,comment);
                }else{
                    // Otherwise, find name and value, seperated by '=' and ';'
                    string name;
                    string val;
                    getline(f,name,'=');
                    getline(f,val,';');
                    name = Trim(name, " \t\n\r");
                    val = Trim(val, " \t\n\r");
                    
                    if( name.size() >0 && val.size() > 0 )
                    {    
						if (name.compare("dx") == 0)
							*dx = stod(val);
						else if (name.compare("dy") == 0)
							*dy = stod(val);
						else if (name.compare("nRows") == 0)
							*nRows = stoi(val);
						else if (name.compare("nCols") == 0)
							*nCols = stoi(val);
						else if (name.compare("model") == 0)
							*model = val;
						else if (name.compare("input_mode") == 0)
							*input_mode = !val.compare("usb_cam");
						else if (name.compare("input_name") == 0)
							input_name = val;
                    }
                }
            }
        }
        f.close();
    }else{
        cerr << "Unable to open '" << filename << "' for configuration data" << endl;
    }
}


