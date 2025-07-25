AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libffarawobjects.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libffarawobjects_la_LIBADD = \
  -lphool

ROOTDICTS = \
  CaloPacket_Dict.cc \
  CaloPacketv1_Dict.cc \
  CaloPacketContainer_Dict.cc \
  CaloPacketContainerv1_Dict.cc \
  Gl1Packet_Dict.cc \
  Gl1Packetv1_Dict.cc \
  Gl1Packetv2_Dict.cc \
  Gl1Packetv3_Dict.cc \
  Gl1RawHit_Dict.cc \
  Gl1RawHitv1_Dict.cc \
  Gl1RawHitv2_Dict.cc \
  InttRawHit_Dict.cc \
  InttRawHitContainer_Dict.cc \
  InttRawHitContainerv1_Dict.cc \
  InttRawHitContainerv2_Dict.cc \
  InttRawHitv1_Dict.cc \
  InttRawHitv2_Dict.cc \
  LL1Packet_Dict.cc \
  LL1Packetv1_Dict.cc \
  LL1PacketContainer_Dict.cc \
  LL1PacketContainerv1_Dict.cc \
  MicromegasRawHit_Dict.cc \
  MicromegasRawHitContainer_Dict.cc \
  MicromegasRawHitContainerv1_Dict.cc \
  MicromegasRawHitContainerv2_Dict.cc \
  MicromegasRawHitContainerv3_Dict.cc \
  MicromegasRawHitv1_Dict.cc \
  MicromegasRawHitv2_Dict.cc \
  MicromegasRawHitv3_Dict.cc \
  MvtxRawHit_Dict.cc \
  MvtxRawHitv1_Dict.cc \
  MvtxRawHitContainer_Dict.cc \
  MvtxRawHitContainerv1_Dict.cc \
  MvtxFeeIdInfo_Dict.cc \
  MvtxFeeIdInfov1_Dict.cc \
  MvtxRawEvtHeader_Dict.cc \
  MvtxRawEvtHeaderv1_Dict.cc \
  MvtxRawEvtHeaderv2_Dict.cc \
  OfflinePacket_Dict.cc \
  OfflinePacketv1_Dict.cc \
  TpcDiode_Dict.cc \
  TpcDiodeContainer_Dict.cc \
  TpcDiodeContainerv1_Dict.cc \
  TpcDiodev1_Dict.cc \
  TpcRawHit_Dict.cc \
  TpcRawHitContainer_Dict.cc \
  TpcRawHitContainerv1_Dict.cc \
  TpcRawHitContainerv2_Dict.cc \
  TpcRawHitContainerv3_Dict.cc \
  TpcRawHitv1_Dict.cc \
  TpcRawHitv2_Dict.cc \
  TpcRawHitv3_Dict.cc

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

pkginclude_HEADERS = \
  CaloPacket.h \
  CaloPacketv1.h \
  CaloPacketContainer.h \
  CaloPacketContainerv1.h \
  Gl1Packet.h \
  Gl1Packetv1.h \
  Gl1Packetv2.h \
  Gl1Packetv3.h \
  Gl1RawHit.h \
  Gl1RawHitv1.h \
  Gl1RawHitv2.h \
  InttRawHit.h \
  InttRawHitContainer.h \
  InttRawHitContainerv1.h \
  InttRawHitContainerv2.h \
  InttRawHitv1.h \
  InttRawHitv2.h \
  LL1Packet.h \
  LL1Packetv1.h \
  LL1PacketContainer.h \
  LL1PacketContainerv1.h \
  MicromegasRawHit.h \
  MicromegasRawHitContainer.h \
  MicromegasRawHitContainerv1.h \
  MicromegasRawHitContainerv2.h \
  MicromegasRawHitContainerv3.h \
  MicromegasRawHitv1.h \
  MicromegasRawHitv2.h \
  MicromegasRawHitv3.h \
  MvtxRawHit.h \
  MvtxRawHitv1.h \
  MvtxRawHitContainer.h \
  MvtxRawHitContainerv1.h \
  MvtxFeeIdInfo.h \
  MvtxFeeIdInfov1.h \
  MvtxRawEvtHeader.h \
  MvtxRawEvtHeaderv1.h \
  MvtxRawEvtHeaderv2.h \
  OfflinePacket.h \
  OfflinePacketv1.h \
  TpcDiode.h \
  TpcDiodeContainer.h \
  TpcDiodeContainerv1.h \
  TpcDiodev1.h \
  TpcRawHit.h \
  TpcRawHitContainer.h \
  TpcRawHitContainerv1.h \
  TpcRawHitContainerv2.h \
  TpcRawHitContainerv3.h \
  TpcRawHitv1.h \
  TpcRawHitv2.h \
  TpcRawHitv3.h

libffarawobjects_la_SOURCES = \
  $(ROOTDICTS) \
  CaloPacketv1.cc \
  CaloPacketContainerv1.cc \
  Gl1Packet.cc \
  Gl1Packetv1.cc \
  Gl1Packetv2.cc \
  Gl1Packetv3.cc \
  Gl1RawHit.cc \
  Gl1RawHitv1.cc \
  Gl1RawHitv2.cc \
  InttRawHitContainerv1.cc \
  InttRawHitContainerv2.cc \
  InttRawHitv1.cc \
  InttRawHitv2.cc \
  LL1Packetv1.cc \
  LL1PacketContainerv1.cc \
  MicromegasRawHitContainerv1.cc \
  MicromegasRawHitContainerv2.cc \
  MicromegasRawHitContainerv3.cc \
  MicromegasRawHitv1.cc \
  MicromegasRawHitv2.cc \
  MicromegasRawHitv3.cc \
  MvtxRawHitv1.cc \
  MvtxRawHitContainerv1.cc \
  MvtxFeeIdInfov1.cc \
  MvtxRawEvtHeaderv1.cc \
  MvtxRawEvtHeaderv2.cc \
  OfflinePacket.cc \
  OfflinePacketv1.cc \
  TpcDiodeContainerv1.cc \
  TpcDiodev1.cc \
  TpcRawHitContainerv1.cc \
  TpcRawHitContainerv2.cc \
  TpcRawHitContainerv3.cc \
  TpcRawHitv1.cc \
  TpcRawHitv2.cc \
  TpcRawHitv3.cc

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libffarawobjects.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

clean-local:
	rm -f *Dict.* *.pcm $(BUILT_SOURCES)
