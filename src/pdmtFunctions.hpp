/*****************************************************************************

         This file is a part of PDMT (Parallel Dual Meshing Tool)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.com
     Date     : 22/02/2022
     Comment  : The program collects functions needed by the PDMT mesher

     -------------------------------------------------------------------

     PDMT a parallel  dual meshing tool uses   finite  element framework
     to convert a triangular / tetrahedral mesh into a  polyhedral  mesh.
     PDMT is distributed  in  the  hope that it  will be useful, HOWEVER
     WITHOUT ANY WARRANTY; or without  even  implied warranty of FITNESS
     FOR A PARTICULAR PURPOSE.

*******************************************************************************/

#include <iostream>
#include <stack>
#include <stdlib.h>

using namespace std;

//-------------------------------------------------------------------------
// PdmtGetMeshInfo :
//     This function takes as inputs i) 2D mesh pTh, ii) vector isBorderNode
//     and ii) a  vector meshInfo. At the end of this function
//          \forall i=[0:Th.nv]
//           isBorderNode(i) = 0 (if i is not on border)
//           isBorderNode(i) > 0 (if i is on border)
//     and
//          meshInfo contains the following info
//            meshInfo(0) = Total # nodes
//            meshInfo(1) = Total # nodes inside
//            meshInfo(2) = Total # nodes on boundary
//            meshInfo(2) = Total # triangles
//            meshInfo(2) = Total # edges
//-------------------------------------------------------------------------
int PdmtGetMeshInfo(const Fem2D::Mesh* const &pTh, KN< long > *const &isBorderNode, KN< long > *const &meshInfo) {

#ifdef DEBUG
  cout << "\n"
          "--------------------------------------\n"
          " PdmtGetBorderInfo function called  \n"
          "--------------------------------------\n";
#endif

  const Mesh &Th = *pTh;

  isBorderNode->resize(Th.nv);                  // resize isBorderNode
  meshInfo->resize(5);

  long nbBorderNodes = 0;

  for(int i = 0; i < Th.nv; ++i){                // isBorderNode = 0
    *(isBorderNode[0]+i) = 0;
//    cout << " Th.nv.lab " << Th(i).lab << endl;
    }

  for(int k = 0; k < Th.neb; ++k)               // isBorderNode > 0
   for(int i = 0; i < 2; ++i){
//     cout << " Th.be(k).lab " << Th.be(k).lab << endl;
     *(isBorderNode[0]+Th(Th.be(k)[i])) += 1;
     }

  for(int i = 0; i < Th.nv; ++i)                // isBorderNode = 0
    if(*(isBorderNode[0]+i))
      nbBorderNodes = nbBorderNodes+1;

  *(meshInfo[0]+0) = Th.nv                 ; // meshInfo(0) = Total # nodes
  *(meshInfo[0]+1) = Th.nv - nbBorderNodes ; // meshInfo(1) = Total # nodes inside
  *(meshInfo[0]+2) = nbBorderNodes         ; // meshInfo(2) = Total # nodes on boundary
  *(meshInfo[0]+3) = Th.nt                 ; // meshInfo(3) = Total # triangles
  *(meshInfo[0]+4) = Th.neb                ; // meshInfo(4) = Total # edges

#ifdef DEBUG
  cout << "\n"
          " PdmtGetBorderInfo stats:\n"
          "  meshInfo[0] - # nodes             = "  << *(meshInfo[0]+0)  << "\n"
          "  meshInfo[1] - # nodes inside      = "  << *(meshInfo[0]+1)  << "\n"
          "  meshInfo[2] - # nodes on boundary = "  << *(meshInfo[0]+2)  << "\n"
          "  meshInfo[3] - # triangles         = "  << *(meshInfo[0]+3)  << "\n"
          "  meshInfo[4] - # edges             = "  << *(meshInfo[0]+4)  << "\n";

  cout << "--------------------------------------\n" << endl;
#endif

  return 0;
}

//-------------------------------------------------------------------------
// PdmtFillSearchTableEdges :
//           This function takes 2D mesh Th and two vectors headVertexBorder
//           and nextVertexBorder as input. At the end these two vectors can
//           be used to search and navigate through number of edges for each
//           mesh point i;
//-------------------------------------------------------------------------
int PdmtFillSearchTableEdges(const Fem2D::Mesh* const &pTh, KN< long > *const &headVertexBorder,  KN< long > *const &nextVertexBorder) {

  const Mesh &Th = *pTh;

#ifdef DEBUG
  cout << "\n"
          "--------------------------------------\n"
          " PdmtFillSearchTableEdges function called  \n"
          "--------------------------------------\n";

  cout << "  Th.nv  - mesh vertices -    "  << Th.nv  << "\n"
          "  Th.nbe - mesh edges    -    "  << Th.neb << "\n";

  cout << "--------------------------------------\n" << endl;
#endif


  headVertexBorder->resize(Th.nv);              // resize headVertexBorder
  nextVertexBorder->resize(Th.neb*2);           // resize nextVertexBorder

  for(int i = 0; i < Th.nv; ++i)                // set headVertexBorder = -1
    *(headVertexBorder[0]+i) = -1;

  for(int k = 0; k < Th.neb; ++k){
    for(int i = 0; i < 2; ++i){
      int v = Th(Th.be(k)[i]);                                        // current vertex number
      *(nextVertexBorder[0]+(2*k+i))  =  *(headVertexBorder[0]+v);    // next vertex
      *(headVertexBorder[0]+v)     = 2*k+i;                           // update head vertex
    }
  }

  return 0;
}


//-------------------------------------------------------------------------
// PdmtFillSearchTableTriangles :
//           This function takes 2D mesh Th and two vectors  headVertex  and
//           and  nextVertex as input. At  the  end these two vectors can be
//           used to search  &  navigate through number of triangles for each
//           mesh point i;
//-------------------------------------------------------------------------
int PdmtFillSearchTableTriangles(const Fem2D::Mesh* const &pTh, KN< long > *const &headVertex,  KN< long > *const &nextVertex) {

  const Mesh &Th = *pTh;

#ifdef DEBUG
  cout << "\n"
          "--------------------------------------\n"
          " PdmtFillSearchTableTriangles function called  \n"
          "--------------------------------------\n";

  cout << "  Th.nv  - mesh vertices -    "  << Th.nv  << "\n"
          "  Th.nbe - mesh edges    -    "  << Th.neb << "\n"
          "  Th.nt  - mesh triangles-    "  << Th.nt  << "\n";

  cout << "--------------------------------------\n" << endl;
#endif


  headVertex->resize(Th.nv);              // resize headVertex
  nextVertex->resize(Th.nt*3);            // resize nextVertex

  for(int i = 0; i < Th.nv; ++i)          // set headVertex = -1
    *(headVertex[0]+i) = -1;

  for(int k = 0; k < Th.nt; ++k){
    for(int i = 0; i < 3; ++i){
      int v = Th(Th[k][i]);                                     // current vertex number
      *(nextVertex[0]+(3*k+i))  =  *(headVertex[0]+v);          // next vertex
      *(headVertex[0]+v)        = 3*k+i;                        // update head vertex
    }
  }

  return 0;
}


double PdmtMaxinTwoP1(KN<double> *const & f, KN<double> *const & f1)
{
  long int nn = f->N();

  for(long int i=0; i<nn; i++)
    *(f[0]+i)=max(*(f1[0]+i),*(f[0]+i));

  return 0.0;
}

//-------------------------------------------------------------------------
// PdmtParseGmshPhysicalTags :
//           Parses named or integer physical line groups from a Gmsh .msh file
//           for viscous layer boundary group selection.
//-------------------------------------------------------------------------
inline std::string pdmtTrimString(const std::string &value) {
  std::string::size_type begin = 0;
  while (begin < value.size() && std::isspace(static_cast<unsigned char>(value[begin])))
    ++begin;
  std::string::size_type end = value.size();
  while (end > begin && std::isspace(static_cast<unsigned char>(value[end - 1])))
    --end;
  return value.substr(begin, end - begin);
}

inline std::vector<std::string> pdmtSplitCsvString(const std::string &csv) {
  std::vector<std::string> result;
  std::string::size_type begin = 0;
  while (begin <= csv.size()) {
    const std::string::size_type comma = csv.find(',', begin);
    const std::string token = pdmtTrimString(csv.substr(
        begin, comma == std::string::npos ? std::string::npos : comma - begin));
    if (!token.empty())
      result.push_back(token);
    if (comma == std::string::npos)
      break;
    begin = comma + 1;
  }
  return result;
}

int PdmtParseGmshPhysicalTags(string* const &pMeshFile,
                              string* const &pRequestedGroups,
                              KN< long > *const &pTags) {
  const std::string meshFileName = *pMeshFile;
  const std::string requestedCsv = *pRequestedGroups;
  const std::vector<std::string> requested = pdmtSplitCsvString(requestedCsv);

  std::ifstream input(meshFileName.c_str());
  std::map<std::string, long> nameToTag;
  std::set<long> allLineTags;

  if (input) {
    std::string line;
    while (std::getline(input, line)) {
      std::string section = pdmtTrimString(line);
      if (section == "$PhysicalNames") {
        std::string countLine;
        if (std::getline(input, countLine)) {
          const long count = std::strtol(countLine.c_str(), 0, 10);
          for (long i = 0; i < count; ++i) {
            std::string pLine;
            if (!std::getline(input, pLine)) break;
            std::istringstream values(pLine);
            long dimension = -1;
            long tag = -1;
            values >> dimension >> tag;
            const std::string::size_type firstQuote = pLine.find('"');
            const std::string::size_type lastQuote = pLine.rfind('"');
            if (dimension == 1 && firstQuote != std::string::npos && lastQuote > firstQuote) {
              std::string pName = pLine.substr(firstQuote + 1, lastQuote - firstQuote - 1);
              nameToTag[pName] = tag;
              allLineTags.insert(tag);
            }
          }
        }
      } else if (section == "$Elements") {
        std::string countLine;
        if (std::getline(input, countLine)) {
          const long count = std::strtol(countLine.c_str(), 0, 10);
          for (long i = 0; i < count; ++i) {
            std::string eLine;
            if (!std::getline(input, eLine)) break;
            std::istringstream values(eLine);
            long id = 0, elementType = 0, numberOfTags = 0;
            values >> id >> elementType >> numberOfTags;
            long physicalTag = 0;
            for (long t = 0; t < numberOfTags; ++t) {
              long val = 0;
              values >> val;
              if (t == 0) physicalTag = val;
            }
            if (elementType == 1 || elementType == 8 || elementType == 26 || elementType == 27 || elementType == 28) {
              if (physicalTag > 0)
                allLineTags.insert(physicalTag);
            }
          }
        }
      }
    }
  }

  std::set<long> matchedTags;
  for (size_t i = 0; i < requested.size(); ++i) {
    const std::string &req = requested[i];
    if (req == "ALL") {
      matchedTags.insert(allLineTags.begin(), allLineTags.end());
    } else if (nameToTag.find(req) != nameToTag.end()) {
      matchedTags.insert(nameToTag[req]);
    } else {
      char *endPtr = 0;
      long val = std::strtol(req.c_str(), &endPtr, 10);
      if (endPtr != req.c_str() && *endPtr == '\0') {
        matchedTags.insert(val);
      }
    }
  }

  pTags->resize(matchedTags.size());
  long idx = 0;
  for (std::set<long>::const_iterator it = matchedTags.begin(); it != matchedTags.end(); ++it) {
    (*pTags)[idx++] = *it;
  }

  return 0;
}

