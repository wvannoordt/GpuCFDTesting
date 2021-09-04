#ifndef OUTPUT_H
#define OUTPUT_H

#include "FlowField.h"
#include <string>
#include "Format.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <cstring>

static inline void Base64ByteConversionStream(std::ostream& strm, char* bytes, size_t size)
{
    //This is a very slow implentation for now.
    const char  table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    char term = '=';
    char b0, b1, b2, c0, c1, c2, c3;
    size_t numWriteTotal = (size - size%3)/3;
    for (size_t pos = 0; pos < numWriteTotal; pos++)
    {
        b0 = bytes[3*pos];
        b1 = bytes[3*pos+1];
        b2 = bytes[3*pos+2];
        c0 = ((b0 & 0b11111100) >> 2);
        c1 = ((b0 & 0b00000011) << 4) | ((b1 & 0b11110000) >> 4);
        c2 = ((b1 & 0b00001111) << 2) | ((b2 & 0b11000000) >> 6);
        c3 = ((b2 & 0b00111111));
        strm << table[c0] << table[c1] << table[c2] << table[c3];
    }
    int numLeft = size%3;
    if (numLeft==0) return;
    char last3[3] = {0};
    char last4[4] = {0};
    int num2write[3] = {0, 2, 3};
    int numTerms[3] = {0, 2, 1};
    for (int i = 0; i < numLeft; i++)
    {
        last3[i] = bytes[size-(size%3)+i];
    }
    b0 = last3[0];
    b1 = last3[1];
    b2 = last3[2];
    c0 = ((b0 & 0b11111100) >> 2);
    c1 = ((b0 & 0b00000011) << 4) | ((b1 & 0b11110000) >> 4);
    c2 = ((b1 & 0b00001111) << 2) | ((b2 & 0b11000000) >> 6);
    c3 = ((b2 & 0b00111111));
    last4[0] = c0;
    last4[1] = c1;
    last4[2] = c2;
    last4[3] = c3;
    for (int i = 0; i < num2write[numLeft]; i++)
    {
        strm << table[last4[i]];
    }
    for (int i = 0; i < numTerms[numLeft]; i++)
    {
        strm << term;
    }
}

static inline bool CreateDirectory(std::string directoryPath)
{
    int check = mkdir(directoryPath.c_str(), 0755);
    if (check != 0)
    {
        auto er = errno;
        if (er==EEXIST) return true;
        std::stringstream ss;
        ss << std::strerror(errno);
        CallError(strformat("System error: {}", ss.str()));
        return false;
    }
    return true;
}

static inline bool MachineIsBigEndian()
{
    int num = 1;
    return (! ( *(char *)&num == 1 ));
}

static inline std::string ZFill(int val, int znum)
{
    std::string output = std::to_string(val);
    while (output.length()<znum) output = "0" + output;
    return output;
}

static inline std::string spaces(int i)
{
    std::string output = "";
    for (int ii = 0; ii < i; ii++) output += " ";
    return output;
}

static void Output(const FlowField& flow, std::string directory, std::string fileTitle)
{
    std::string blockDirectory = strformat("{}/{}_b", directory, fileTitle);
    std::string mainFilename   = strformat("{}/{}.vtm", directory, fileTitle);
    std::string blockTemplateFileNameRelative = strformat("{}_b/", fileTitle) + "b{}.vtr";
    CreateDirectory(blockDirectory);
    
    std::string blockTemplateFileName = blockDirectory + "/b{}.vtr";
    
    
    bool is3D = flow.is3D;
    for (int lb = 0; lb < flow.numBlocks; lb++)
    {
        print("Block", lb, "of", flow.numBlocks-1);
        int nCellsi = flow.imax - 2*flow.nguard;
        int nCellsj = flow.jmax - 2*flow.nguard;
        int nCellsk = flow.kmax - 2*flow.nguard;
        if (!is3D) nCellsk = 1;
        int nGuardi = flow.nguard;
        int nGuardj = flow.nguard;
        int nGuardk = flow.nguard;
        if (!is3D) nGuardk = 0;
        int nTotali = nCellsi + 2*nGuardi;
        int nTotalj = nCellsj + 2*nGuardj;
        int nTotalk = nCellsk + 2*nGuardk;
        
        
        
        std::string filename = strformat(blockTemplateFileName, ZFill(lb, 7));
        auto box = flow.GetBox(lb);
        double ghostBnds[6];
        ghostBnds[4] = 0.0;
        ghostBnds[5] = 0.0;
        for (int i = 0; i < (is3D?3:2); i++)
        {
            ghostBnds[2*i] = box.bounds[2*i] - box.dx[i]*flow.nguard;
            ghostBnds[2*i+1] = box.bounds[2*i+1] + box.dx[i]*flow.nguard;
        }
        std::string varsString = "";
        for (int l = 0; l < flow.varNames.size(); l++)
        {
            varsString += ((l==0)?"":",");
            varsString += flow.varNames[l];
        }
        
        std::ofstream myfile;
        myfile.open(filename.c_str());
        std::string sp20 = spaces(20);
        const char* csp20 = sp20.c_str();
        std::string endianness = MachineIsBigEndian()?"BigEndian":"LittleEndian";
        myfile << "<?xml version=\"1.0\"?>\n";
        myfile << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"" << endianness << "\" header_type=\"UInt32\">" << std::endl;
        myfile << spaces(4) << strformat("<RectilinearGrid WholeExtent=\"0 {} 0 {} 0 {}\">", nTotali, nTotalj, nTotalk) << std::endl;
        myfile << spaces(8) << "<FieldData>" << std::endl;
        myfile << spaces(12) << "<DataArray type=\"Int32\" Name=\"avtRealDims\" NumberOfTuples=\"6\" format=\"ascii\">" << std::endl;
        myfile << spaces(16) << strformat("{} {} {} {} {} {}", nGuardi, nGuardj+nCellsi, nGuardj, nGuardj+nCellsj, nGuardk, is3D?(nGuardk+nCellsk):0) << std::endl;
        myfile << spaces(12) << "</DataArray>" << std::endl;
        myfile << spaces(12) << "<DataArray type=\"Float64\" Name=\"avtOriginalBounds\" NumberOfTuples=\"6\" format=\"ascii\">" << std::endl;
        myfile << spaces(16) << strformat("{} {} {} {} {} {}", box.bounds[0], box.bounds[1], box.bounds[2], box.bounds[3], box.bounds[4], box.bounds[5]) << std::endl;
        myfile << spaces(12) << "</DataArray>" << std::endl;
        myfile << spaces(8) << "</FieldData>" << std::endl;
        myfile << spaces(8) << strformat("<Piece Extent=\"0 {} 0 {} 0 {}\">", nTotali, nTotalj, nTotalk) << std::endl;
        
        myfile << spaces(12) << "<CellData Scalars=\"" << varsString << "\">" << std::endl;
        
        auto array = flow.UnloadBlock(lb);
        for (int var = 0; var < flow.numVars; var++)
        {
            myfile << spaces(16) << "<DataArray type=\"Float64\" Name=\"" << flow.varNames[var] << "\" format=\"binary\">" << std::endl;
            ////https://mathema.tician.de/what-they-dont-tell-you-about-vtk-xml-binary-formats/
            unsigned int ss = array.dims[0]*array.dims[1]*array.dims[2]*sizeof(double);
            size_t offset = ss;
            Base64ByteConversionStream(myfile, (char*)(&ss), sizeof(unsigned int));
            Base64ByteConversionStream(myfile, (char*)(array.data)+var*offset, offset);
            myfile << "\n" << spaces(16) << "</DataArray>" << std::endl;
        }
        
        myfile << spaces(16) << "<DataArray type=\"UInt8\" Name=\"avtGhostZones\" format=\"ascii\">" << std::endl;
        auto isGhost = [&](int i, int j, int k) -> bool
        {
            return
                (i<flow.nguard)||
                (i>=array.dims[0]-flow.nguard)||
                (j<flow.nguard)||
                (j>=array.dims[1]-flow.nguard)||
                (k<flow.nguard)||
                (k>=array.dims[2]-flow.nguard);
        };
        for (int k = 0; k < flow.kmax; k++)
        {
            for (int j = 0; j < flow.jmax; j++)
            {
                for (int i = 0; i < flow.imax; i++)
                {
                    myfile << csp20 << (isGhost(i, j, k)?16:0) << "\n";
                }
            }
        }
        myfile << spaces(16) << "</DataArray>" << std::endl;
        myfile << spaces(12) << "</CellData>" << std::endl;
        myfile << spaces(12) << "<Coordinates>" << std::endl;
        myfile << spaces(16) << strformat("<DataArray type=\"Float64\" format=\"ascii\" RangeMin=\"{}\" RangeMax=\"{}\">", ghostBnds[0], ghostBnds[1]) << std::endl;
        for (int i = -nGuardi; i <=nCellsi+nGuardi; i++)
        {
            myfile << csp20 << box.bounds[0] + i*box.dx[0] << "\n";
        }
        myfile << spaces(16) << "</DataArray>" << std::endl;
        myfile << spaces(16) << strformat("<DataArray type=\"Float64\" format=\"ascii\" RangeMin=\"{}\" RangeMax=\"{}\">", ghostBnds[2], ghostBnds[3]) << std::endl;
        for (int j = -nGuardj; j <=nCellsj+nGuardj; j++)
        {
            myfile << csp20 << box.bounds[2] + j*box.dx[1] << "\n";
        }
        myfile << spaces(16) << "</DataArray>" << std::endl;
        myfile << spaces(16) << strformat("<DataArray type=\"Float64\" format=\"ascii\" RangeMin=\"{}\" RangeMax=\"{}\">", ghostBnds[4], ghostBnds[4]) << std::endl;
        for (int k = -nGuardk; k <=nCellsk+nGuardk; k++)
        {
            myfile << csp20 << (is3D?(box.bounds[4] + k*box.dx[2]):0.0) << "\n";
        }
        myfile << spaces(16) << "</DataArray>" << std::endl;
        myfile << spaces(12) << "</Coordinates>" << std::endl;
        myfile << spaces(8) << "</Piece>" << std::endl;
        myfile << spaces(4) << "</RectilinearGrid>" << std::endl;
        myfile << "</VTKFile>" << std::endl;
        
        myfile.close();
    }
    
    std::string filename = mainFilename;
    std::ofstream myMfile;
    myMfile.open(filename.c_str());
    myMfile << "<?xml version=\"1.0\"?>\n";
    myMfile << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\">" << std::endl;
    myMfile << spaces(4) << "<vtkMultiBlockDataSet>" << std::endl;
    myMfile << spaces(8) << "<Block index =\"0\">" << std::endl;
    for (int b = 0; b < flow.numBlocks; b++)
    {
        std::string blockFileName = strformat(blockTemplateFileNameRelative, ZFill(b, 7));
        myMfile << spaces(12) << strformat("<DataSet index=\"{}\" file=\"{}\"/>", b, blockFileName) << std::endl;
    }
    myMfile << spaces(8) << "</Block>" << std::endl;
    myMfile << spaces(4) << "</vtkMultiBlockDataSet>" << std::endl;
    myMfile << "</VTKFile>" << std::endl;
    myMfile.close();
}

#endif