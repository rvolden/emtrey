package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	sc "strconv"
	str "strings"
	"time"
	uni "unicode"
)

// command line arguments
var in = flag.String("i", "", "Input SAM file.")
var mm = flag.Bool("m", false, "Use if your SAM file is from minimap2.")

func check(e error) {
	if e != nil {
		panic(e)
	}
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
	const maxBuffer = 1024 * 1024
	buf := make([]byte, 0, maxBuffer)
	scanner.Buffer(buf, maxBuffer)

	// keep hash of chr:size
	chroms := make(map[string]string)
	for scanner.Scan() {
		line := scanner.Text()
		if line[0] == '@' { // SAM header lines
			if line[1] != 'S' {
				continue
			}
			sLine := str.Split(line, "\t")
			// putting chr sizes into map
			chroms[sLine[1][3:]] = sLine[2][3:]
		} else { // alignment lines
			sLine := str.Split(line, "\t")
			name, qSize := sLine[0], chroms[sLine[2]]
			if sLine[2] == "*" {
				continue
			}
			bitFlag, _ := sc.Atoi(sLine[1])
			var strand string
			if (uint(bitFlag) & 0x10) != 0 {
				strand = "-"
			} else {
				strand = "+"
			}

			tStart, _ := sc.Atoi(sLine[3])
			tStart--

			splitC := parseCIGAR(sLine[5])
			var qStart, qEnd, end, M, I, D, N, S, H, EQ, X, nI, nD int
			var blockSizes = []string{}
			var qStarts = []int{}
			var tStarts = []int{tStart}
			for j, entry := range splitC {
				last := len(entry) - 1
				// adjust start of the read according the clipping
				if j == 0 && str.Contains("SH", string(entry[last])) {
					qStart, _ = sc.Atoi(entry[:last])
					qStarts = append(qStarts, qStart)
				} else if j == 0 {
					qStarts = append(qStarts, 0)
				}
				// also adjust the end
				if j == len(splitC)-1 && entry[last] == 'S' {
					qEnd, _ = sc.Atoi(entry[:last])
				}
				switch entry[last] {
				case 'M': // matches and mismatches
					m, _ := sc.Atoi(entry[:last])
					M += m
					blockSizes = append(blockSizes, entry[:last])
					qStarts = append(qStarts, m+qStarts[len(qStarts)-1])
					tStarts = append(tStarts, m+tStarts[len(tStarts)-1])
				case 'I': // insertions
					i, _ := sc.Atoi(entry[:last])
					I += i
					nI++
					qStarts[len(qStarts)-1] += i
				case 'D': // deletions
					d, _ := sc.Atoi(entry[:last])
					D += d
					nD++
					tStarts[len(tStarts)-1] += d
				case 'N': // unmapped bases
					n, _ := sc.Atoi(entry[:last])
					N += n
					tStarts[len(tStarts)-1] += n
				case 'S': // soft clipped bases
					s, _ := sc.Atoi(entry[:last])
					S += s
				case 'H': // hard clipped bases
					h, _ := sc.Atoi(entry[:last])
					H += h
				case '=': // matches
					eq, _ := sc.Atoi(entry[:last])
					EQ += eq
				case 'X': // mismatches
					x, _ := sc.Atoi(entry[:last])
					X += x
				}
			}

			ID := I + D                    // indels
			sLen := M + I + S + H + EQ + X // calculated sequence lenght
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
				for _, col := range sLine[9:] {
					if str.Contains(col, "NM:i:") {
						NM, _ = sc.Atoi(col[5:])
					} else if str.Contains(col, "nn:i:") {
						ambig, _ = sc.Atoi(col[5:])
					} else if str.Contains(col, "ts:A:") {
						newStrand := col[5:]
						if newStrand == "-" && strand == "+" {
							strand = "-"
						} else if newStrand == "-" && strand == "-" {
							strand = "+"
						}
					}
				}
				mismatch = NM - ID - ambig
				if mismatch < 0 {
					mismatch = 0
				}
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
