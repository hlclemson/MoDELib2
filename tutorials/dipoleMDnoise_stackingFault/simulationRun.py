import re, os
import subprocess
import itertools
import shutil
import numpy as np
from dataclasses import dataclass
from collections import defaultdict
from typing import Any


@dataclass
class dislocationDynamicsRun:
    structure: object
    testRange: dict

    def exploreAllParams(self) -> None:
        modelibPath = self.structure.configFile["mainMoDELibDirectory"]
        inputFilePath = f"{modelibPath}/tutorials/DislocationDynamics/periodicDomains/uniformLoadController/inputFiles/"
        forceFilePath = f"{modelibPath}/tutorials/DislocationDynamics/periodicDomains/uniformLoadController/F/"
        outputPath = self.structure.configFile["dataOutPutDirectory"]
        paramToCouple = self.structure.configFile["paramtersToCouple"]
        testTimeStep = self.structure.configFile["testTimeSteps"]
        # totalTimeStep = self.structure.configFile['totalTimeSteps']
        microStruct = self.structure.configFile["microstructureFileToUse"]
        workingSimPath = f"{modelibPath}/tutorials/DislocationDynamics/periodicDomains/uniformLoadController/"
        microStructLibPath = (
            f"{modelibPath}/tutorials/DislocationDynamics/MicrostructureLibrary"
        )
        materialLibPath = (
            f"{modelibPath}/tutorials/DislocationDynamics/MaterialsLibrary"
        )
        noiseLibPath = f"{modelibPath}/tutorials/DislocationDynamics/NoiseLibrary"
        externalLoadMode = self.structure.configFile["loadType"]
        slipSystemType = self.structure.configFile["slipSystemType"]
        # stressRange = self.structure.configFile['stressTestRange']
        stackingFaultNoise = self.structure.configFile["stackingFaultNoiseMode"]
        solidSolutionNoise = self.structure.configFile["solidSolutionNoiseMode"]
        whiteNoise = self.structure.configFile["whiteNoiseMode"]
        extBoundaryCondition = self.structure.configFile["extBCsToTest"]
        extBCsign = self.structure.configFile["extBCsign"]
        glidePlasticStrain = float(self.structure.configFile["glidePlasticStrain"])

        # remove the dictionary emlements that are empty
        keysToRemove = [
            key for key, values in self.testRange.items() if not values
        ]  # keys to remove
        for key in keysToRemove:  # remove the keys that has empty values
            del self.testRange[key]

        tempKeys = []
        tempArray = []
        # get the keys and values from the testRange dictionary
        for key, value in self.testRange.items():
            tempKeys.append(key)
            tempArray.append(value)

        # create a list of dictionaries for each simulation run
        paramDictList = []
        combinations = list(itertools.product(*tempArray))
        for comb in combinations:
            templateRunDict = {}
            for i, element in enumerate(comb):
                templateRunDict[tempKeys[i]] = element
            paramDictList.append(templateRunDict)

        # clean up the old data in the simulation working directory
        if os.path.exists(f"{workingSimPath}/F"):
            os.system(f"rm -rf {workingSimPath}/F")
        if os.path.exists(f"{workingSimPath}/evl"):
            os.system(f"rm -rf {workingSimPath}/evl")
        # clean up the old data in the data output directory
        if os.path.exists(f"{outputPath}/F"):
            os.system(f"rm -rf {outputPath}/F")
        if os.path.exists(f"{outputPath}/evl"):
            os.system(f"rm -rf {outputPath}/evl")

        # run simulations with the parameters saved on each list
        for parameters in paramDictList:
            print(f"currently running simulation with parameters... : {parameters}")
            # create stress range to test
            # extBCs = np.array(extBoundaryCondition)
            currentSeed = parameters["seed"]
            # find the external boundary conditions assgined to the current seed
            for seedAssgin, extbc in extBoundaryCondition.items():
                if seedAssgin == f"seed_{currentSeed}":
                    #extBC = np.array(extbc)
                    extBC = extbc
                    break
            # extBCs = np.array()
            #if not extBCs.any():
            if not extBC:
                print(
                    f"Boundary condition values for the seed {currentSeed} are not specified."
                )
                continue

            # extract material info
            b_SI = self.readValFromMaterialFile("b_SI", materialLibPath, parameters)
            mu0_SI = self.readValFromMaterialFile("mu0_SI", materialLibPath, parameters)
            rho_SI = self.readValFromMaterialFile("rho_SI", materialLibPath, parameters)
            cs = np.sqrt(mu0_SI / rho_SI)  # shear wave speed
            convertTimeUnit = b_SI / cs  # [sec]
            convertMPaToMu = 1 / (mu0_SI * 10 ** (-6))
            unitConvFactors = {
                "convertTimeUnit": convertTimeUnit,
                "convertMPaToMu": convertMPaToMu,
            }

            # change parameters
            _ = self.changeParameters(
                parameters,
                stackingFaultNoise,
                solidSolutionNoise,
                whiteNoise,
                modelibPath,
                inputFilePath,
                microStructLibPath,
            )
            # if partial is enabled, make the change on the material file
            _ = self.setSlipSystemType(parameters, materialLibPath, slipSystemType)
            # set time step
            _ = self.setTimeStep(testTimeStep, modelibPath, inputFilePath)

            # start searching
            searchStepSize = 100
            searchResolution = 10
            dislocMovingSlowly = False
            dislocationMoving = False
            while not dislocMovingSlowly:
            #for extBC in extBCs:
                print(f"testing boundary condition, extBC = {extBC}")
                sigma = self.modifyExternalLoad(
                    inputFilePath,
                    convertTimeUnit,
                    convertMPaToMu,
                    externalLoadMode,
                    extBC,
                    extBCsign,
                )

                # generate microstructure
                _ = self.generateMicrostructure(modelibPath)

                # run simulation for the number of steps in config.json file
                _ = self.runDislocationDynamics(modelibPath, externalLoadMode)

                # check if the dislocation is moving or not
                # read labels
                with open(f"{workingSimPath}/F/F_labels.txt", "r") as label:
                    fLabels = label.read()
                # remove empty element and store it as a list
                fLabels = np.array([x for x in fLabels.split("\n") if x], dtype=str)

                # last 33 element indexes
                lastCols = np.arange(-33, 0)
                # open F file
                fData = defaultdict(list)
                with open(f"{workingSimPath}/F/F_0.txt", "r") as f:
                    for line in f:
                        line = [float(x) for x in line.split(" ") if x and x != "\n"]
                        # parse the first 14 elements
                        for i in range(len(fLabels[:14])):
                            fData[str(fLabels[i])].append(line[i])
                        # parse the last 33 elements
                        for j in lastCols:
                            fData[str(fLabels[j])].append(line[j])

                xAxisData = np.array(fData["runID"])
                yAxisData = np.array(fData["betaP_13"])
                # absolute plastic strain
                yAxisData = np.abs(yAxisData - yAxisData[0])

                dydx = np.gradient(yAxisData, xAxisData)
                # return the mean p-strain rate of the last 5 recorded steps
                #dydx_mean = np.mean(dydx[-5:])
                # return the mean p-strain rate of the last 3 recorded steps
                dydx_mean = np.mean(dydx[-3:])

                if dydx_mean > glidePlasticStrain:
                    # while loop stop condition
                    #if searchStepSize < 10 and dislocationMoving:
                    if searchStepSize < searchResolution:
                        dislocMovingSlowly = True
                    else:
                        dislocationMoving = True
                        searchStepSize -= searchStepSize//2
                        extBC -= searchStepSize
                else:
                    if searchStepSize < searchResolution and dislocationMoving:
                        dislocMovingSlowly = True
                    elif dislocationMoving:
                        searchStepSize -= searchStepSize//2
                        extBC += searchStepSize
                    else:
                        extBC += searchStepSize

                # stop the run if the boundary condition becomes less than zero
                if extBC < 0:
                    break

                # save the data
                # _ = self.copyDataToOutputDir(parameters, outputPath, workingSimPath, microStructLibPath, noiseLibPath, parameters['alloy'], microStruct, externalLoadMode, stackingFaultNoise, solidSolutionNoise, sigma/convertMPaToMu)
                _ = self.copyDataToOutputDir(
                    parameters,
                    outputPath,
                    workingSimPath,
                    microStructLibPath,
                    noiseLibPath,
                    parameters["alloy"],
                    microStruct,
                    externalLoadMode,
                    stackingFaultNoise,
                    solidSolutionNoise,
                    sigma,
                    **unitConvFactors,
                )

                # remove old file
                if os.path.exists(f"{workingSimPath}/F"):
                    os.system(f"rm -rf {workingSimPath}/F")
                if os.path.exists(f"{workingSimPath}/evl"):
                    os.system(f"rm -rf {workingSimPath}/evl")

    def modifyExternalLoad(
        self,
        inputFilePath: str,
        convertTimeUnit: float,
        convertMPaToMu: float,
        externalLoadMode: str,
        extBC: float,
        extBCsign: str,
    ) -> float:
        match externalLoadMode:
            case "ExternalStress0":
                # read the uniformExternalLoadController.txt file
                with open(
                    f"{inputFilePath}/uniformExternalLoadController.txt", "r"
                ) as file:
                    text = file.read()
                # set the first stress as the lower bound
                sigma = extBC * convertMPaToMu

                # set the rest to zero
                pattern = r"ExternalStrain0.=((?:.|\s)*?);"
                replace = f"ExternalStrain0 = 0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0;"
                # replace the pattern with the new value
                text = re.sub(pattern, replace, text)

                # set the rest to zero
                pattern = r"ExternalStrainRate.=((?:.|\s)*?);"
                replace = f"ExternalStrainRate = 0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0;"
                # replace the pattern with the new value
                text = re.sub(pattern, replace, text)

                # set the rest to zero
                pattern = r"ExternalStressRate.=((?:.|\s)*?);"
                replace = f"ExternalStressRate = 0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0;"
                text = re.sub(pattern, replace, text)

                # reset the stiffness
                pattern = r"MachineStiffnessRatio.=((?:.|\s)*?);"
                replace = f"MachineStiffnessRatio =0.0 0.0 0.0 0.0 0.0 0.0;"
                # replace the pattern with the new value
                text = re.sub(pattern, replace, text)

                # set the initial stress
                pattern = r"ExternalStress0.=((?:.|\s)*?);"
                replace = f"ExternalStress0 = 0.0 0.0 {extBCsign}{sigma}\n0.0 0.0 0.0\n{extBCsign}{sigma} 0.0 0.0;"
                # replace the pattern with the new value
                text = re.sub(pattern, replace, text)

                # overwrite the uniformExternalLoadController.txt file with the new stress
                with open(
                    f"{inputFilePath}/uniformExternalLoadController.txt", "w"
                ) as file:
                    file.write(text)
            case "ExternalStressRate":
                # read the uniformExternalLoadController.txt file
                with open(
                    f"{inputFilePath}/uniformExternalLoadController.txt", "r"
                ) as file:
                    text = file.read()

                # set the rest to zero
                pattern = r"ExternalStrain0.=((?:.|\s)*?);"
                replace = f"ExternalStrain0 = 0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0;"
                # replace the pattern with the new value
                text = re.sub(pattern, replace, text)

                # set the rest to zero
                pattern = r"ExternalStrainRate.=((?:.|\s)*?);"
                replace = f"ExternalStrainRate = 0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0;"
                # replace the pattern with the new value
                text = re.sub(pattern, replace, text)

                # set the rest to zero
                pattern = r"ExternalStress0.=((?:.|\s)*?);"
                replace = f"ExternalStress0 = 0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0;"
                text = re.sub(pattern, replace, text)

                # reset the stiffness
                pattern = r"MachineStiffnessRatio.=((?:.|\s)*?);"
                replace = f"MachineStiffnessRatio =0.0 0.0 0.0 0.0 0.0 0.0;"
                # replace the pattern with the new value
                text = re.sub(pattern, replace, text)

                pattern = r"ExternalStressRate.=((?:.|\s)*?);"
                replace = f"ExternalStressRate = 0.0 0.0 {extBCsign}{sigma}\n0.0 0.0 0.0\n{extBCsign}{sigma} 0.0 0.0;"
                # replace the pattern with the new value
                text = re.sub(pattern, replace, text)

                # overwrite the uniformExternalLoadController.txt file with the new stress
                with open(
                    f"{inputFilePath}/uniformExternalLoadController.txt", "w"
                ) as file:
                    file.write(text)
            case "ExternalStrain0":
                # read the uniformExternalLoadController.txt file
                with open(
                    f"{inputFilePath}/uniformExternalLoadController.txt", "r"
                ) as file:
                    text = file.read()

                # set the rest to zero
                pattern = r"ExternalStrainRate.=((?:.|\s)*?);"
                replace = f"ExternalStrainRate = 0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0;"
                text = re.sub(pattern, replace, text)

                # set the rest to zero
                pattern = r"ExternalStress0.=((?:.|\s)*?);"
                replace = f"ExternalStress0 = 0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0;"
                text = re.sub(pattern, replace, text)

                # set the rest to zero
                pattern = r"ExternalStressRate.=((?:.|\s)*?);"
                replace = f"ExternalStressRate = 0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0;"
                text = re.sub(pattern, replace, text)

                # reset the stiffness
                pattern = r"MachineStiffnessRatio.=((?:.|\s)*?);"
                replace = f"MachineStiffnessRatio =0.0 0.0 0.0 0.0 0.0 0.0;"
                # replace the pattern with the new value
                text = re.sub(pattern, replace, text)

                pattern = r"ExternalStrain0.=((?:.|\s)*?);"
                replace = f"ExternalStrain0 = 0.0 0.0 {extBCsign}{sigma}\n0.0 0.0 0.0\n{extBCsign}{sigma} 0.0 0.0;"
                # replace the pattern with the new value
                text = re.sub(pattern, replace, text)

                # overwrite the uniformExternalLoadController.txt file with the new stress
                with open(
                    f"{inputFilePath}/uniformExternalLoadController.txt", "w"
                ) as file:
                    file.write(text)
            case "ExternalStrainRate":
                # read the uniformExternalLoadController.txt file
                with open(
                    f"{inputFilePath}/uniformExternalLoadController.txt", "r"
                ) as file:
                    text = file.read()
                # rate[1/sec]*ddTime[sec/(b/cs)] = ddRate[1/(b/cs)]
                sigma = extBC * convertTimeUnit

                # set the rest to zero
                pattern = r"ExternalStress0.=((?:.|\s)*?);"
                replace = f"ExternalStress0 = 0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0;"
                text = re.sub(pattern, replace, text)

                # set the rest to zero
                pattern = r"ExternalStrain0.=((?:.|\s)*?);"
                replace = f"ExternalStrain0 = 0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0;"
                # replace the pattern with the new value
                text = re.sub(pattern, replace, text)

                # set the rest to zero
                pattern = r"ExternalStressRate.=((?:.|\s)*?);"
                replace = f"ExternalStressRate = 0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0;"
                text = re.sub(pattern, replace, text)

                pattern = r"ExternalStrainRate.=((?:.|\s)*?);"
                replace = f"ExternalStrainRate = 0.0 0.0 {extBCsign}{sigma}\n0.0 0.0 0.0\n{extBCsign}{sigma} 0.0 0.0;"
                # replace the pattern with the new value
                text = re.sub(pattern, replace, text)

                # set the stiffness
                pattern = r"MachineStiffnessRatio.=((?:.|\s)*?);"
                replace = f"MachineStiffnessRatio =0.0 0.0 0.0 0.0 0.0 1e20;"
                # replace the pattern with the new value
                text = re.sub(pattern, replace, text)

                # overwrite the uniformExternalLoadController.txt file with the new stress
                with open(
                    f"{inputFilePath}/uniformExternalLoadController.txt", "w"
                ) as file:
                    file.write(text)
        return sigma

    def saveData(self, parameters: dict, outputPath: str, sigmaMPa: float) -> None:
        # data output name
        dataOutputName = "CRSStestResult.txt"
        header, data = [], []
        for key, value in parameters.items():
            # if the value is string with ' ' delimiter, split with delim, then strip newline (\n) if there is any
            value = (
                [x.strip("\n") for x in value.split(" ")]
                if type(value) == str
                else value
            )
            # if the value is list, join the string with ',' delimiter
            value = ",".join([str(x) for x in value]) if type(value) == list else value
            match key:
                case "temperature":
                    header.append(f"Temp")
                    data.append(f"{value}")
                case "alloy":
                    header.append(f"Alloy")
                    data.append(f"{value}")
                case "boxSize":
                    header.append(f"BoxSize")
                    data.append(f"{value}")
                case "alphaLineTension":
                    header.append(f"LineTension")
                    data.append(f"{value}")
                case "periodicDipoleSlipSystemIDs":
                    header.append(f"sIDs")
                    data.append(f"{value}")
                case "periodicDipoleExitFaceIDs":
                    header.append(f"exIDs")
                    data.append(f"{value}")
                case "periodicDipoleNodes":
                    header.append(f"NodeNum")
                    data.append(f"{value}")
                case "periodicDipolePoints":
                    header.append(f"dipolePoints")
                    data.append(f"{value}")
                case "periodicDipoleHeights":
                    header.append(f"dipoleHeights")
                    data.append(f"{value}")
                case "periodicDipoleGlideSteps":
                    header.append(f"dipoleGSteps")
                    data.append(f"{value}")
                case "dxMax":
                    header.append(f"dxMax")
                    data.append(f"{value}")
        # add CRSS header
        header.append(f"CRSS")
        # add CRSS
        data.append(f"{sigmaMPa}")

        # join the list of strings as a single string
        header = " ".join(header)
        data = " ".join(data)
        # if file does not exsit, write header
        if not os.path.exists(f"{outputPath}/{dataOutputName}"):
            with open(f"{outputPath}/{dataOutputName}", "w") as output:
                output.write(f"{header}\n")
                output.write(f"{data}\n")
        else:
            with open(f"{outputPath}/{dataOutputName}", "a") as output:
                output.write(f"{data}\n")

    def readValFromMaterialFile(
        self, parameter: str, libPath: str, parameters: dict
    ) -> float:
        with open(f"{libPath}/{parameters['alloy']}.txt", "r") as mFile:
            for line in mFile:
                # strip down comments from the data
                line = line.split(";")[0]
                if line.startswith(f"{parameter}"):
                    # Split the line by '=' and take the second part, then remove leading/trailing whitespace
                    value = float(line.split("=")[1].strip())
                    break
        return value

    def copyDataToOutputDir(
        self,
        parameters: dict,
        outputPath: str,
        workingSimPath: str,
        microStructLibPath: str,
        noiseLibPath: str,
        alloy: str,
        microStructFile: str,
        externalLoadMode: str,
        stackingFaultNoise: int,
        solidSolutionNoise: int,
        sigma: float,
        **kwargs,
    ) -> None:
        # create output directory if there isn't one
        if not os.path.exists(f"{outputPath}"):
            shutil.os.makedirs(f"{outputPath}")

        convertTimeUnit = kwargs["convertTimeUnit"]
        convertMPaToMu = kwargs["convertMPaToMu"]
        # create directory name based on the parameters
        name = []
        match externalLoadMode:
            case "ExternalStress0":
                name.append(f"Str{sigma/convertMPaToMu:.0f}")
                print(f"Str{sigma/convertMPaToMu:.0f}")
            case "ExternalStressRate":
                name.append(f"StrR{sigma*convertTimeUnit:.0f}")
            case "ExternalStrain0":
                name.append(f"Strn{sigma:.0f}")
            case "ExternalStrainRate":
                print(f"StrnR{sigma/convertTimeUnit:.0e}")
                name.append(f"StrnR{sigma/convertTimeUnit:.0e}")
        # name.append(f'Str{int(sigmaMPa)}')
        for key, value in parameters.items():
            # if the value is string with ' ' delimiter, split with delim, then strip newline (\n) if there is any
            value = (
                [x.strip("\n") for x in value.split(" ")]
                if type(value) == str
                else value
            )
            # if the value is list, join the string with '' delimiter
            value = "".join([str(x) for x in value]) if type(value) == list else value
            match key:
                case "temperature":
                    name.append(f"T{value}")
                case "alloy":
                    name.append(f"{value}")
                case "boxSize":
                    name.append(f"BS{value}")
                case "alphaLineTension":
                    name.append(f"LT{value}")
                case "periodicDipoleSlipSystemIDs":
                    name.append(f"sID{value}")
                case "periodicDipoleExitFaceIDs":
                    name.append(f"exID{value}")
                case "periodicDipoleNodes":
                    name.append(f"N{value}")
                case "periodicDipolePoints":
                    name.append(f"DP{value}")
                case "periodicDipoleHeights":
                    name.append(f"DH{value}")
                case "periodicDipoleGlideSteps":
                    name.append(f"DGS{value}")
                case "dxMax":
                    name.append(f"dMx{value}")
                case "seed":
                    name.append(f"S{value}")
        # join the list of strings as a single string
        folderName = "".join(name)

        # Remove the old generated data if there is previously generated data
        if os.path.exists(f"{outputPath}/{folderName}"):
            os.system(f"rm -rf {outputPath}/{folderName}")

        # copy the data to the output directory
        shutil.copytree(
            f"{workingSimPath}/", f"{outputPath}/{folderName}", dirs_exist_ok=True
        )
        # copy microStructureFile
        shutil.copy(
            f"{microStructLibPath}/{microStructFile}", f"{outputPath}/{folderName}"
        )

        if solidSolutionNoise:
            # copy noise files
            with open(f"{workingSimPath}/inputFiles/polycrystal.txt", "r") as file:
                polyTxt = file.read()
            noiseNamesInPoly = re.findall("solidSolutionNoiseFile.*", polyTxt)
            ssNoiseToCopy = []
            for noiseName in noiseNamesInPoly:
                # print(f'noiseName = {noiseName.split('/')[-1].strip(';')}')
                ssNoiseToCopy.append(noiseName.split("/")[-1].strip(";"))
            # copy sampled solid solution noise file
            # shutil.copy(f'{noiseLibPath}/noise_xz.vtk', f'{outputPath}/{folderName}/noise_xz_SS_{alloy}.vtk')
            # shutil.copy(f'{noiseLibPath}/noise_yz.vtk', f'{outputPath}/{folderName}/noise_yz_SS_{alloy}.vtk')
            for ssNoiseFile in ssNoiseToCopy:
                shutil.copy(
                    f"{noiseLibPath}/{ssNoiseFile}", f"{outputPath}/{folderName}/"
                )
        if stackingFaultNoise:
            # copy noise files
            with open(f"{workingSimPath}/inputFiles/polycrystal.txt", "r") as file:
                polyTxt = file.read()
            noiseNamesInPoly = re.findall("stackingFaultNoiseFile.*", polyTxt)
            sfNoiseToCopy = []
            for noiseName in noiseNamesInPoly:
                # print(f'noiseName = {noiseName.split('/')[-1].strip(';')}')
                sfNoiseToCopy.append(noiseName.split("/")[-1].strip(";"))
            # copy sampled stacking fault noise file
            # shutil.copy(f'{noiseLibPath}/noise_{alloy}.vtk', f'{outputPath}/{folderName}/noiseSF_{alloy}.vtk')
            for sfNoiseFile in sfNoiseToCopy:
                shutil.copy(
                    f"{noiseLibPath}/{sfNoiseFile}", f"{outputPath}/{folderName}/"
                )

        # clean up the data in the working directory
        if os.path.exists(f"{workingSimPath}/F"):
            os.system(f"rm -rf {workingSimPath}/F")
        if os.path.exists(f"{workingSimPath}/evl"):
            os.system(f"rm -rf {workingSimPath}/evl")

    def generateMicrostructure(self, modelibPath: str) -> None:
        # modelibPath = self.structure.configFile['mainMoDELibDirectory']
        binaryFile = (
            f"{modelibPath}/tools/MicrostructureGenerator/build/microstructureGenerator"
        )
        workingSimPath = f"{modelibPath}/tutorials/DislocationDynamics/periodicDomains/uniformLoadController/"
        necessaryFolders = ["evl", "F"]
        for folder in necessaryFolders:
            # Check if the directory does not exist, create it
            if not os.path.exists(f"{workingSimPath}/{folder}"):
                os.makedirs(f"{workingSimPath}/{folder}")
        # catch the binary runtime error
        try:
            # Execute the binary
            result = subprocess.run(
                [f"{binaryFile}", f"{workingSimPath}"],
                check=True,  # Raises a CalledProcessError on non-zero exit status
                stdout=subprocess.PIPE,  # Capture standard output
                stderr=subprocess.PIPE,  # Capture standard error
            )
            # output = result.stdout.decode('utf-8')
            error = result.stderr.decode("utf-8")
            # print("Output:", output)
            # print("Error:", error)
        except subprocess.CalledProcessError as e:
            exit(e.stderr.decode("utf-8"))

    def runDislocationDynamics(self, modelibPath: str, externalLoadMode: str) -> None:
        # print(f'currently running simulation with parameters... : {parameters}')
        binaryFile = f"{modelibPath}/tools/DDomp/build/DDomp"
        workingSimPath = f"{modelibPath}/tutorials/DislocationDynamics/periodicDomains/uniformLoadController/"
        inputFilePath = f"{modelibPath}/tutorials/DislocationDynamics/periodicDomains/uniformLoadController/inputFiles/"

        max_attempts = 3  # Maximum number of attempts to restart the binary
        attempt = 1
        while attempt <= max_attempts:
            try:
                # Execute the binary
                result = subprocess.run(
                    [f"{binaryFile}", f"{workingSimPath}"],
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                error = result.stderr.decode("utf-8")
                break  # Exit the loop if the execution is successful
            except subprocess.CalledProcessError as e:
                error = e.stderr.decode("utf-8")
                print(f"Attempt {attempt} failed. Error: {error}")
                attempt += 1

        # don't stop the entire simulation
        if attempt > max_attempts:
            print(f"Execution failed after {max_attempts} attempts.")
        #    exit(error)

    def calcMeanOfDotBetaP(self, forceFilePath: str, convertTimeUnit: float) -> float:
        try:
            # read the data
            time = 1
            s13bPIndex = 5
            tstep = np.loadtxt(f"{forceFilePath}/F_0.txt", usecols=time)
            data = np.loadtxt(f"{forceFilePath}/F_0.txt", usecols=s13bPIndex)

            # Calculate the index to start from (50% of the array length)
            startIndex = int(len(data) * 0.5)

            # Slice the array to get the remaining 80%
            trimmedData = data[startIndex:]
            trimmedTstep = tstep[startIndex:]

            # calculate the first derivative of the remaining 80% of the time-plasticStrain graph
            diff = np.gradient(
                trimmedData, trimmedTstep * convertTimeUnit
            )  # rate [1/s]
            muDiff = np.mean(diff)  # mean of the plastic strain rate [1/s]
            return muDiff
        except FileNotFoundError:
            exit(
                f"File not found: {forceFilePath}/F_0.txt\n Check if DDomp is executed properly"
            )
        except ValueError as ve:
            exit(f"ValueError occurred: {str(ve)}")
        except Exception as e:
            exit(f"An unexpected error occurred: {str(e)}")

    def copyReferenceInputFiles(self, modelibPath: str) -> None:
        # Define the source and destination directories
        sourceDir = f"./ReferenceInputFiles/"
        # Copy a directory and its contents, overwriting the destination if it already exists
        shutil.copytree(
            f"{sourceDir}",
            f"{modelibPath}/tutorials/DislocationDynamics/",
            dirs_exist_ok=True,
        )

    def setTimeStep(self, timeStep: int, modelibPath: str, inputFilePath: str) -> None:
        ddFile = f"{inputFilePath}/DD.txt"
        pattern = f"Nsteps=.*"
        replace = f"Nsteps={timeStep};"
        with open(ddFile, "r") as file:
            text = file.read()
        # replace the pattern with the new value
        text = re.sub(pattern, replace, text)
        # overwrite the original data file
        with open(ddFile, "w") as file:
            file.write(text)

    def modifyTXTfile(self, filePath: str, pattern: str, replace: str) -> None:
        with open(filePath, "r") as file:
            text = file.read()
        # replace the pattern with the new value
        text = re.sub(pattern, replace, text)
        # overwrite the original data file
        with open(filePath, "w") as file:
            file.write(text)

    def modifyGenInputFilePy(
        self, filePath: str, pattern: str, replace: str, inputFilePath: str
    ) -> None:
        with open(filePath, "r") as file:
            text = file.read()
        # replace the pattern with the new value
        text = re.sub(pattern, replace, text)
        # overwrite the original data file
        with open(filePath, "w") as file:
            file.write(text)
        # change the current working directory to inputFiles
        runtimeDir = os.getcwd()
        os.chdir(f"{inputFilePath}")
        # using exec() to run another Python script
        with open(f"./generateInputFiles.py", "r") as file:
            exec(file.read())
        # change the current working directory back to the original
        os.chdir(runtimeDir)

    def setSlipSystemType(
        self, paramDictionary: dict, materialLibPath: str, slipSystemType: str
    ) -> None:
        for key, value in paramDictionary.items():
            match key:
                case "alloy":
                    pattern = f"enabledSlipSystems.*"
                    replace = f"enabledSlipSystems={slipSystemType};"
                    filePath = f"{materialLibPath}/{value}.txt"
        with open(filePath, "r") as file:
            text = file.read()
        # replace the pattern with the new value
        if re.search(pattern, text):
            text = re.sub(pattern, replace, text)
            # overwrite the original data file
            with open(filePath, "w") as file:
                file.write(text)
        # if there is no enabledSlipSystem declaration, add it to the file
        else:
            with open(filePath, "w") as file:
                file.write(text)
                file.write(f"{replace}")

    def changeParameters(
        self,
        paramDictionary: dict,
        stackingFaultNoise: int,
        solidSolutionNoise: int,
        whiteNoise: int,
        modelibPath: str,
        inputFilePath: str,
        microStructLibPath: str,
    ) -> None:
        # change solidSolutionNoiseMode
        pattern = f"pf.solidSolutionNoiseMode.*"
        replace = f"pf.solidSolutionNoiseMode={solidSolutionNoise};"
        filePath = f"{inputFilePath}/generateInputFiles.py"
        self.modifyGenInputFilePy(filePath, pattern, replace, inputFilePath)
        # change stackingFaultNoiseMode
        pattern = f"pf.stackingFaultNoiseMode.*"
        replace = f"pf.stackingFaultNoiseMode={stackingFaultNoise}; # 0=no noise"
        filePath = f"{inputFilePath}/generateInputFiles.py"
        self.modifyGenInputFilePy(filePath, pattern, replace, inputFilePath)
        # change parameters
        for key, value in paramDictionary.items():
            match key:
                case "temperature":
                    pattern = f"pf.absoluteTemperature.*"
                    replace = f"pf.absoluteTemperature={value};"
                    filePath = f"{inputFilePath}/generateInputFiles.py"
                    self.modifyGenInputFilePy(filePath, pattern, replace, inputFilePath)
                case "alloy":
                    pattern = f"pf=PolyCrystalFile.*"
                    replace = (
                        f'pf=PolyCrystalFile("../../../MaterialsLibrary/{value}.txt");'
                    )
                    filePath = f"{inputFilePath}/generateInputFiles.py"
                    self.modifyGenInputFilePy(filePath, pattern, replace, inputFilePath)
                    if whiteNoise:
                        pattern = f"pf.stackingFaultCorrelationFile.*"
                        replace = f"pf.stackingFaultCorrelationFile='../../../NoiseLibrary/{value}_Cx_R100_ISF_w.vtk';"
                        filePath = f"{inputFilePath}/generateInputFiles.py"
                        self.modifyGenInputFilePy(
                            filePath, pattern, replace, inputFilePath
                        )
                        pattern = f"pf.stackingFaultNoiseFile.*"
                        replace = f"pf.stackingFaultNoiseFile='../../../NoiseLibrary/sf_noise_{value}_w.vtk';"
                        filePath = f"{inputFilePath}/generateInputFiles.py"
                        self.modifyGenInputFilePy(
                            filePath, pattern, replace, inputFilePath
                        )
                    else:
                        pattern = f"pf.stackingFaultCorrelationFile.*"
                        replace = f"pf.stackingFaultCorrelationFile='../../../NoiseLibrary/{value}_Cx_R100_ISF.vtk';"
                        filePath = f"{inputFilePath}/generateInputFiles.py"
                        self.modifyGenInputFilePy(
                            filePath, pattern, replace, inputFilePath
                        )
                        pattern = f"pf.stackingFaultNoiseFile.*"
                        replace = f"pf.stackingFaultNoiseFile='../../../NoiseLibrary/sf_noise_{value}.vtk';"
                        filePath = f"{inputFilePath}/generateInputFiles.py"
                        self.modifyGenInputFilePy(
                            filePath, pattern, replace, inputFilePath
                        )
                # case 'lineTension':
                #    pattern = f'alphaLineTension=.*'
                #    replace = f'alphaLineTension={value};'
                #    filePath = f'{inputFilePath}/DD.txt'
                #    self.modifyTXTfile(filePath, pattern, replace)
                case "boxSize":
                    pattern = f"pf.boxScaling.*"
                    replace = f"pf.boxScaling=np.array({value});"
                    filePath = f"{inputFilePath}/generateInputFiles.py"
                    self.modifyGenInputFilePy(filePath, pattern, replace, inputFilePath)
                case "seed":
                    pattern = f"pf.seed.*"
                    replace = f"pf.seed={value};"
                    filePath = f"{inputFilePath}/generateInputFiles.py"
                    self.modifyGenInputFilePy(filePath, pattern, replace, inputFilePath)
                case "solidSolutionGridSpacing_SI":
                    pattern = f"pf.solidSolutionGridSpacing_SI.*"
                    replace = f"pf.solidSolutionGridSpacing_SI=np.array({value});"
                    filePath = f"{inputFilePath}/generateInputFiles.py"
                    self.modifyGenInputFilePy(filePath, pattern, replace, inputFilePath)
                case "solidSolutionGridSize":
                    pattern = f"pf.solidSolutionGridSize.*"
                    replace = f"pf.solidSolutionGridSize=np.array({value});"
                    filePath = f"{inputFilePath}/generateInputFiles.py"
                    self.modifyGenInputFilePy(filePath, pattern, replace, inputFilePath)
                case "stackingFaultGridSize":
                    pattern = f"pf.stackingFaultGridSize.*"
                    replace = f"pf.stackingFaultGridSize=np.array({value});"
                    filePath = f"{inputFilePath}/generateInputFiles.py"
                    self.modifyGenInputFilePy(filePath, pattern, replace, inputFilePath)
                case "coreSize":
                    pattern = f"pf.a_cai_A.*"
                    replace = f"pf.a_cai_A={float(value)};"
                    filePath = f"{inputFilePath}/generateInputFiles.py"
                    self.modifyGenInputFilePy(filePath, pattern, replace, inputFilePath)
                case "periodicDipoleSlipSystemIDs":
                    pattern = f"periodicDipoleSlipSystemIDs.*"
                    replace = f"periodicDipoleSlipSystemIDs={value};"
                    filePath = f"{microStructLibPath}/periodicDipole.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "periodicDipoleExitFaceIDs":
                    pattern = f"periodicDipoleExitFaceIDs.*"
                    replace = f"periodicDipoleExitFaceIDs={value};"
                    filePath = f"{microStructLibPath}/periodicDipole.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "periodicDipoleNodes":
                    pattern = f"periodicDipoleNodes.*"
                    replace = f"periodicDipoleNodes={value};"
                    filePath = f"{microStructLibPath}/periodicDipole.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "periodicDipolePoints":
                    pattern = r"periodicDipolePoints=((?:.|\s)*?);"
                    replace = f"periodicDipolePoints={value};"
                    filePath = f"{microStructLibPath}/periodicDipole.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "periodicDipoleHeights":
                    pattern = f"periodicDipoleHeights.*"
                    replace = f"periodicDipoleHeights={value};"
                    filePath = f"{microStructLibPath}/periodicDipole.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "periodicDipoleGlideSteps":
                    pattern = f"periodicDipoleGlideSteps.*"
                    replace = f"periodicDipoleGlideSteps={value};"
                    filePath = f"{microStructLibPath}/periodicDipole.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "minimizationSteps":
                    pattern = f"relaxSteps.*"
                    replace = f"relaxSteps={value};"
                    filePath = f"{inputFilePath}/uniformExternalLoadController.txt"
                    self.modifyTXTfile(filePath, pattern, replace)

        # change parameters in DD.txt files, it needs to be done in this order separately because the python script overwrites the txt files
        for key, value in paramDictionary.items():
            match key:
                case "alphaLineTension":
                    pattern = f"alphaLineTension=.*"
                    replace = f"alphaLineTension={value};"
                    filePath = f"{inputFilePath}/DD.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "coreSize":
                    pattern = f"coreSize.*"
                    replace = f"coreSize={value};"
                    filePath = f"{inputFilePath}/DD.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "dxMax":
                    pattern = f"dxMax.*"
                    replace = f"dxMax={value};"
                    filePath = f"{inputFilePath}/DD.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "periodicImages":
                    pattern = f"periodicImageSize.*"
                    replace = f"periodicImageSize={value};"
                    filePath = f"{inputFilePath}/DD.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "outputFrequency":
                    pattern = f"outputFrequency.*"
                    replace = f"outputFrequency={value};"
                    filePath = f"{inputFilePath}/DD.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "EwaldLengthFactor":
                    pattern = f"EwaldLengthFactor.*"
                    replace = f"EwaldLengthFactor={value};"
                    filePath = f"{inputFilePath}/DD.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "outputSlipSystemPlasticDistortion":
                    pattern = f"outputSlipSystemPlasticDistortion.*"
                    replace = f"outputSlipSystemPlasticDistortion={value};"
                    filePath = f"{inputFilePath}/DD.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "maxJunctionIterations":
                    pattern = f"maxJunctionIterations.*"
                    replace = f"maxJunctionIterations={value};"
                    filePath = f"{inputFilePath}/DD.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "timeIntegrationMethod":
                    pattern = f"timeIntegrationMethod.*"
                    replace = f"timeIntegrationMethod={value};"
                    filePath = f"{inputFilePath}/DD.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "timeStep":
                    pattern = f"timeStep.*"
                    replace = f"timeStep={value};"
                    filePath = f"{inputFilePath}/DD.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "quadPerLength":
                    pattern = f"quadPerLength.*"
                    replace = f"quadPerLength={value};"
                    filePath = f"{inputFilePath}/DD.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
                case "computeDDinteractions":
                    pattern = f"computeDDinteractions.*"
                    replace = f"computeDDinteractions={value};"
                    filePath = f"{inputFilePath}/DD.txt"
                    self.modifyTXTfile(filePath, pattern, replace)
