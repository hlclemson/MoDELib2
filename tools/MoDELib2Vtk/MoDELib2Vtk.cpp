/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <MoDELib2Vtk.h>

using namespace model;

int main(int argc, char** argv)
{
#ifdef _MODEL_PYBIND11_ // COMPILED WITH PYBIND11
    pybind11::scoped_interpreter guard{};
#endif

    const std::string folderName(argc>1? std::string(argv[1]) : "./");

    // create an object
    DD2OvitoVtk dd2ov(folderName);

    return 0;
}

