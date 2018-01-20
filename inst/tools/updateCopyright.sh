#!/bin/sh

cd ./src;
echo "Update src header files";
find . -name "*.h" -exec sed -i "s/Copyright 2007-2017 The OpenMx/Copyright 2007-2018 The OpenMx/g" '{}' \;
echo "Update src cpp files";
find . -name "*.cpp" -exec sed -i "s/Copyright 2007-2017 The OpenMx/Copyright 2007-2018 The OpenMx/g" '{}' \;
echo "Update R files";
cd ../R;
find . -name "*.R" -exec sed -i "s/Copyright 2007-2017 The OpenMx/Copyright 2007-2018 The OpenMx/g" '{}' \;
echo "Update man pages";
cd ../man
find . -name "*.Rd" -exec sed -i "s/Copyright 2007-2017 The OpenMx/Copyright 2007-2018 The OpenMx/g" '{}' \;
echo "Update demo R files";
cd ../demo
find . -name "*.R" -exec sed -i "s/Copyright 2007-2017 The OpenMx/Copyright 2007-2018 The OpenMx/g" '{}' \;
echo "Update models/enormous files";
cd ../inst/models/enormous;
find . -name "*.R" -exec sed -i "s/Copyright 2007-2017 The OpenMx/Copyright 2007-2018 The OpenMx/g" '{}' \;
echo "Update models/failing files";
cd ../failing;
find . -name "*.R" -exec sed -i "s/Copyright 2007-2017 The OpenMx/Copyright 2007-2018 The OpenMx/g" '{}' \;
echo "Update models/nightly files";
cd ../nightly;
find . -name "*.R" -exec sed -i "s/Copyright 2007-2017 The OpenMx/Copyright 2007-2018 The OpenMx/g" '{}' \;
echo "Update models/passing files";
cd ../passing;
find . -name "*.R" -exec sed -i "s/Copyright 2007-2017 The OpenMx/Copyright 2007-2018 The OpenMx/g" '{}' \;
cd ../../..;
