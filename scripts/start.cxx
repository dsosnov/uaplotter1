{
gROOT->ProcessLine(".include ../UATree/UADataFormat/interface");
gROOT->ProcessLine(".include ../TOTEMdataFormat/interface");
gROOT->ProcessLine(".include ../interface");
gROOT->ProcessLine(".L ../UATree/UADataFormat/lib/libUADataFormat.so");
gROOT->ProcessLine(".L ../TOTEMdataFormat/lib/libTOTEMdataFormat.so");
gROOT->ProcessLine(".L lib/libuaplotter.so");
}
