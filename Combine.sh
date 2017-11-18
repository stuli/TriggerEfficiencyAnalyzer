#!/bin/bash -f


indir="/eos/cms/store/group/phys_heavyions/jaebeom/OniaTree_PPref_SingleMuTnP/SingleMuonTnP/OniaTree_PPref_SingleMuTnP/171116_031052/0001/"
indir2="/eos/cms/store/group/phys_heavyions/jaebeom/OniaTree_PPref_SingleMuTnP/SingleMuonTnP/OniaTree_PPref_SingleMuTnP/171116_031052/0000/"

ls $indir | grep root > list
grep "root" ./list | awk -v p=$indir '{print "fcha->Add(\""p$1"\");"}' > list2

ls $indir2 | grep root > list0
grep "root" ./list0 | awk -v p=$indir2 '{print "fcha->Add(\""p$1"\");"}' > list3

cat list2 list3 > mergedlist

sed -i '82r ./mergedlist' getEffAllData_1117.cpp




