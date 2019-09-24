BEGIN { 
    FS="=";
    unite[0] = "o";
    unite[1] = "Ko";
    unite[2] = "Mo";
    unite[3] = "Go";
}
    $1 == "mem_heap_B" {
        snapshots[line] = $2; 
    }
    $1 == "mem_heap_extra_B" {
        snapshots[line] += $2; 
    }
    $1 == "mem_stacks_B" {
        snapshots[line] += $2; 
        if (snapshots[line] > maxi) {
            maxi = snapshots[line];
            maxi_line = line
        }
        line = line + 1;
    }

END {
    choix = 0;

    while (maxi > 1000 && choix < 3) {
        maxi = maxi / 1000;
        choix++;
    }
    print(maxi, unite[choix]);
}