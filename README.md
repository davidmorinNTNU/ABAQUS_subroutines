# ABAQUS subroutines

After 15+ years of coding in ABAQUS, I decided to put up a public repository of the subroutines I am using as part of my research. These subroutines (in reality a small part of what is possible in ABAQUS) covers some aspects of material mechanics.

>While I took care of verifying my code (verification) they might not be bug free and you should use them in accordance. Furthermore, the proposed implementation are not necessarily validated (validation) in the sense that the chosen models are merely examples and their performance will be highly case dependent.

In this repository, you will find the following subroutines:

1. UVARM
2. USDFLD / VUSDFLD
3. UHARD / VUHARD
4. UMAT / VUMAT
5. UEXTERNALDB / VEXTERNALDB

In addition to these subroutines, there is an example where a VUHARD subroutine and a VUSDFLD subroutine are combined to implement the Modified Johnson-Cook model.

Each example/folder has the following content:

1. a FORTRAN subroutine,
2. an ABAQUS input file fo both Standard and Explicit (except for UVARM),
3. a python post-processing script to extract results from ABAQUS/viewer,
4. a python script to plot the results,
5. a pdf presentation of the subroutine, the input file structure and the results.

>All analyses have been checked in ABAQUS 2019 with double precision. The Intel compiler ifort version 16.0.8 was employed. While the pdf presentations are (hopefully) useful, the official ABAQUS documentation should always be read first.

To execute the ABAQUS jobs with the FORTRAN subroutines you can use the following command:

`abaqus double job=INPUT_FILE user=USER_SUBROUTINE int`

where: `INPUT_FILE` is the name of the ABAQUS input file and `USER_SUBROUTINE` is the name of the FORTRAN file. For windows users, you might have to change the extension from `*.f` to `*.for`.