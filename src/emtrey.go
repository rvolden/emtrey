package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	sc "strconv"
	str "strings"
	uni "unicode"
	"time"
)

// command line arguments
var in = flag.String("i", "", "Input sam file.")
var mm = flag.Bool("m", false, "Use if your SAM file is from minimap2.")

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func isSet(n, k uint) bool {
	k--
	if n&(1<<k) != 0 {
		return true
	}
	return false
}

func parseCIGAR(cigar string) (splitC []string) {
	// splits CIGAR string into a splice with the number + letter
	tmp := ""
	for _, c := range cigar {
		if uni.IsDigit(c) {
			tmp += string(c)
		} else {
			tmp += string(c)
			splitC = append(splitC, tmp)
			tmp = ""
		}
	}
	return splitC
}

func readSAM(samFile *string) {
	file, err := os.Open(*samFile)
	check(err)
	defer file.Close()

	scanner := bufio.NewScanner(file)

	// keep hash of chr:size
	chroms := make(map[string]string)
	for scanner.Scan() {
		line := scanner.Text()
		if line[0] == '@' {
			if line[1] != 'S' {
				continue
			}
			sLine := str.Split(line, "\t")
			// putting chr sizes into map
			chroms[sLine[1][3:]] = sLine[2][3:]
		} else {
			sLine := str.Split(line, "\t")
			name, qSize := sLine[0], chroms[sLine[2]]

			bitFlag, _ := sc.Atoi(sLine[1])
			var strand string
			if isSet(uint(bitFlag), 16) {
				strand = "-"
			} else {
				strand = "+"
			}

			tStart, _ := sc.Atoi(sLine[3])
			tStart--

			splitC := parseCIGAR(sLine[5])
			var qStart, qEnd, end, M, I, D, N, S, EQ, X, nI, nD int
			var blockSizes = []string{}
			var qStarts = []int{}
			var tStarts = []int{tStart}
			for j, entry := range splitC {
				last := len(entry) - 1
				if j == 0 && str.Contains("SH", string(entry[last])) {
					qStart, _ = sc.Atoi(entry[:last])
					qStarts = append(qStarts, qStart)
				}
				if j == len(splitC)-1 && str.Contains("SH", string(entry[last])) {
					qEnd, _ = sc.Atoi(entry[:last])
				}
				switch entry[last] {
				case 'M':
					m, _ := sc.Atoi(entry[:last])
					M += m
					blockSizes = append(blockSizes, entry[:last])
					qStarts = append(qStarts, m+qStarts[len(qStarts)-1])
					tStarts = append(tStarts, m+tStarts[len(tStarts)-1])
				case 'I':
					i, _ := sc.Atoi(entry[:last])
					I += i
					nI++
					qStarts[len(qStarts)-1] += i
				case 'D':
					d, _ := sc.Atoi(entry[:last])
					D += d
					nD++
					tStarts[len(tStarts)-1] += d
				case 'N':
					n, _ := sc.Atoi(entry[:last])
					N += n
					tStarts[len(tStarts)-1] += n
				case 'S':
					s, _ := sc.Atoi(entry[:last])
					S += s
				case '=':
					eq, _ := sc.Atoi(entry[:last])
					EQ += eq
				case 'X':
					x, _ := sc.Atoi(entry[:last])
					X += x
				}
			}

			ID := I + D
			sLen := M + I + S + EQ + X
			consumeRef := M + D + N + EQ + X
			tEnd := tStart + consumeRef
			if qEnd == 0 {
				end = sLen
			} else {
				end = sLen - qEnd
			}

			if len(qStarts) > 0 {
				qStarts = qStarts[:len(qStarts)-1]
			}
			if len(tStarts) > 0 {
				tStarts = tStarts[:len(tStarts)-1]
			}

			blockCount := len(blockSizes)

			// use minimap's extra flags to determine mismatches
			var NM, ambig, matches, mismatch int
			if *mm {
				for _, col := range sLine {
					if str.Contains(col, "NM:i:") {
						NM, _ = sc.Atoi(col[5:])
					}
					if str.Contains(col, "nn:i:") {
						ambig, _ = sc.Atoi(col[5:])
					}
					if str.Contains(col, "ts:A:") {
						strand = col[5:]
					}
				}
				mismatch = NM - ID - ambig
				matches = M - mismatch
			} else if EQ != 0 {
				matches, mismatch = EQ, X
			} else {
				mismatch, matches = 0, M
			}
			bSize := str.Trim(str.Join(str.Fields(fmt.Sprint(blockSizes)), ","), "[]")
			qSt := str.Trim(str.Join(str.Fields(fmt.Sprint(qStarts)), ","), "[]")
			tSt := str.Trim(str.Join(str.Fields(fmt.Sprint(tStarts)), ","), "[]")

			fmt.Printf("%d\t%d\t0\t%d\t%d\t%d\t", matches, mismatch, N, nI, I)
			fmt.Printf("%d\t%d\t%v\t%v\t%v\t", nD, D, strand, name, sLen)
			fmt.Printf("%d\t%d\t%v\t%v\t%d\t", qStart, end, sLine[2], qSize, tStart)
			fmt.Printf("%d\t%d\t%v,\t%v,\t%v,\n", tEnd, blockCount, bSize, qSt, tSt)
		}
	}
}

func main() {
	start := time.Now()
	flag.Parse()

	if *in == "" {
		log.Fatal("Please specify an input file.\n")
	}
	readSAM(in)

	fmt.Fprintf(os.Stderr, "SAM from mm: %v\n", *mm)

	// stop timer and print elapsed time to stderr
	stop := time.Now()
	elapsed := stop.Sub(start)
	fmt.Fprintf(os.Stderr, "Took %v to run.\n", elapsed)
}
