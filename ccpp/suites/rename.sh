#!/bin/bash
oldname=$1
newname=$2
sed -i.backup "s/$oldname/$newname/g" suite_${oldname}.xml || exit 1
gsed -i "2i <!--\nDESCRIPTION:\nSuite \"$newname\" (formerly known as $oldname)\n-->" suite_${oldname}.xml || exit 2
git mv suite_${oldname}.xml ${newname}.xml
git add ${newname}.xml
echo "    \"${oldname}.xml\": \"${newname}.xml\"," >> alias.json
