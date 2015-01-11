#ifndef __DEFINE__
#define __DEFINE__
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstdlib>
#include <algorithm>
#include "mkl.h"
#include <float.h>
#include <math.h>
#include <iomanip>
#include "arrays.hpp"
#include "structs.hpp"
#include "mkl_lapack.h"
#include "mkl_vsl.h"
#include "mkl_types.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"

class TableHelp
{
public:
	TableHelp()
	:pagewidth_(80)
	{}
	void setpagewidth(int pw) const {pagewidth_ = pw;}
	std::string line() const;
	std::string bigline() const;
	std::string starline() const;
private:
	mutable int pagewidth_;
};
#endif
