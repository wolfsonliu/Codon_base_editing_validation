#! /usr/bin/bash

function splitfq_insertion {
    awk '$6 ~ /(I|S)/ {
        print $0;
    }' $1
}

function splitfq_deletion {
    awk '$6 ~ /D/ {
        print $0;
    }' $1
}

function splitfq_other {
    awk '$6 !~ /(I|S|D)/ {
        print $0;
    }' $1
}

for x in CZZ58_S57_L008 CZZ59_S58_L008 CZZ8_S6_L002; do
    splitfq_insertion ./result/${x}.sorted.sam > ./result/${x}.insertion.sam
    splitfq_deletion ./result/${x}.sorted.sam > ./result/${x}.deletion.sam
    splitfq_other ./result/${x}.sorted.sam > ./result/${x}.other.sam
done
