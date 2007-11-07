all: gen_model_cords gcproj

getpar:
        @cd ../Getpar/getpar/src && ${MAKE}

gen_model_cords: getpar
        @cd GenGrid_py/src && ${MAKE}

gcproj: getpar
        @cd Modelbox/src && ${MAKE}

clean:
        @cd ../Getpar/getpar/src && ${MAKE} clean
        @cd GenGrid_py/src && ${MAKE} clean
        @cd Modelbox/src && ${MAKE} clean
