#include "WriteData.h"



template<class Coord>
void Write_Radiation_vtr(double* E, double* Fx, double* Fy, Grid<Coord>& grid, int nFrames, std::string directory)
{
	int NX = grid.n1;
	int NY = grid.n2;
	int NZ = 1;
	for(unsigned int f=0; f<nFrames; f++)
	{
		std::string Filename = directory + "/xdata_" + FrameNumber(f) + ".vtr";
		std::ofstream file(Filename);

		file << "<?xml version=\"1.0\"?>" << std::endl;
		file << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
		file << "<RectilinearGrid WholeExtent=\""
			 << 0 << " " << NX-1 << " " << 0 << " " << NY-1 << " " << 0  << " " << NZ-1 << "\">" << std::endl;
		file << "\t<Piece Extent=\""
			 << 0 << " " << NX-1 << " " << 0 << " " << NY-1 << " " << 0  << " " << NZ-1 << "\">" << std::endl;
		file << "\t\t<PointData>" << std::endl;

		file << "\t\t\t<DataArray type=\"Float32\" Name=\"Moment 0\" format=\"ascii\">" << std::endl;
		for(int ij=0; ij<NX*NY; ij++)
			file << "\t\t\t\t" << E[ij + f*NX*NY*NZ] << std::endl;
		file << "\t\t\t</DataArray>" << std::endl;
		file << "\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"2\" Name=\"Moment 1\" format=\"ascii\">" << std::endl;
		for(int ij=0; ij<NX*NY; ij ++)
			file << "\t\t\t\t" << Fx[ij + f*NX*NY*NZ] << " " << Fy[ij + f*NX*NY*NZ] << std::endl;
		file << "\t\t\t</DataArray>" << std::endl;

		file << "\t\t</PointData>"  << std::endl;
		file << "\t\t<CellData />"  << std::endl;

		file << "\t\t<Coordinates>" << std::endl;
		file << "\t\t\t<DataArray type=\"Float32\" Name=\"X\" format=\"ascii\">" << std::endl;
		file << "\t\t\t\t";
		for(int i=0; i<NX; i++)
			file << grid.xCoord(i,0) << " ";
		file << std::endl;
		file << "\t\t\t</DataArray>" << std::endl;
		file << "\t\t\t<DataArray type=\"Float32\" Name=\"Y\" format=\"ascii\">" << std::endl;
		file << "\t\t\t\t";
		for(int j=0; j<NY; j++)
			file << grid.yCoord(0,j) << " ";
		file << std::endl;
		file << "\t\t\t</DataArray>" << std::endl;
		file << "\t\t\t<DataArray type=\"Float32\" Name=\"Z\" format=\"ascii\">" << std::endl;
		file << "\t\t\t\t";
		for(int k=0; k<NZ; k++)
			file << k+1 << " ";
		file << std::endl;
		file << "\t\t\t</DataArray>" << std::endl;
		file << "\t\t</Coordinates>" << std::endl;

		file << "\t</Piece>" << std::endl;
		file << "</RectilinearGrid>" << std::endl;
		file << "</VTKFile>" << std::endl;
		file.close();
	}
}



template<class Coord>
void Write_Radiation_json(double* E, double* Fx, double* Fy, Grid<Coord>& grid, int nFrames, std::string name)
{
	std::string Filename = name + ".json";
	std::ofstream file(Filename);
	std::string base_filename = name.substr(name.find_last_of("/\\") + 1);
	// Data overhead:
	file << "{\n";
	file << "\t\"name\": \"" << base_filename << "\",\n";
	file << "\t\"dimensions\": 2,\n";
	file << "\t\"numberOfFrames\": " << nFrames << ",\n";
	file << "\t\"n1\": " << grid.n1 << ",\n";
	file << "\t\"n2\": " << grid.n2 << ",\n";
	file << "\t\"d1\": " << grid.d1 << ",\n";
	file << "\t\"d2\": " << grid.d2 << ",\n";
	file << "\t\"start1\": " << grid.start1 << ",\n";
	file << "\t\"start2\": " << grid.start2 << ",\n";
	file << "\t\"end1\": " << grid.end1 << ",\n";
	file << "\t\"end2\": " << grid.end2 << ",\n";
	file << "\t\"frames\": [\n";

	for(unsigned int f=0; f<nFrames; f++)
	{
		file << "\t\t{\n";

		// Data E:
		file << "\t\t\t\"E\": [";
		for(int ij=0; ij<grid.n12; ij++)
		{
			if(isnan(E[ij + f*grid.n12]))
				file << 0;
			else
				file << E[ij + f*grid.n12];
			if(ij<grid.n12-1)
				file << ",";
		}
		file << "],\n";

		// Data Fx:
		file << "\t\t\t\"Fx\": [";
		for(int ij=0; ij<grid.n12; ij++)
		{
			if(isnan(Fx[ij + f*grid.n12]))
				file << 0;
			else
				file << Fx[ij + f*grid.n12];
			if(ij<grid.n12-1)
				file << ",";
		}
		file << "],\n";

		// Data Fy:
		file << "\t\t\t\"Fy\": [";
		for(int ij=0; ij<grid.n12; ij++)
		{
			if(isnan(Fy[ij + f*grid.n12]))
				file << 0;
			else
				file << Fy[ij + f*grid.n12];
			if(ij<grid.n12-1)
				file << ",";
		}
		file << "]\n";

		file << "\t\t}";
		if(f<nFrames-1)
			file << ",\n";
		else
			file << "\n";
	}
	file << "\t]\n";
	file << "}";
}



template void Write_Radiation_vtr<xy>  (double* E, double* Fx, double* Fy, Grid<xy>&  grid, int nFrames, std::string directory);
template void Write_Radiation_vtr<rph> (double* E, double* Fx, double* Fy, Grid<rph>& grid, int nFrames, std::string directory);
template void Write_Radiation_json<xy> (double* E, double* Fx, double* Fy, Grid<xy>&  grid, int nFrames, std::string name);
template void Write_Radiation_json<rph>(double* E, double* Fx, double* Fy, Grid<rph>& grid, int nFrames, std::string name);