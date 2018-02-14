{
string s_include_dir = string(".include ");
string s_library_dir = string(".L ");
string cmssw_base = gSystem->Getenv("CMSSW_BASE");

gROOT->ProcessLine( (s_include_dir+cmssw_base+string("/src/uaplotter1/interface")).c_str() );
gROOT->ProcessLine( (s_include_dir+cmssw_base+string("/src/UATree/UADataFormat/interface")).c_str() );
gROOT->ProcessLine( (s_include_dir+cmssw_base+string("/src/TOTEMdataFormat/interface")).c_str() );

gROOT->ProcessLine((s_library_dir+cmssw_base+string("/src/UATree/UADataFormat/lib/libUADataFormat.so")).c_str());
gROOT->ProcessLine((s_library_dir+cmssw_base+string("/src/TOTEMdataFormat/lib/libTOTEMdataFormat.so")).c_str());
gROOT->ProcessLine((s_library_dir+cmssw_base+string("/src/uaplotter1/lib/libuaplotter.so")).c_str());
}
