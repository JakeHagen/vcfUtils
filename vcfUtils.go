package main

import (
	"bufio"
	"context"
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/brentp/faidx"
	"github.com/brentp/vcfgo"
	"github.com/google/subcommands"
)

type manipInfo struct {
	fields   string
	prefix   string
	operator string
}

func (*manipInfo) Name() string     { return "manipInfo" }
func (*manipInfo) Synopsis() string { return "Make new info field based off other fields" }
func (*manipInfo) Usage() string {
	return `manipInfo -operator [max,min,mean] -prefix prefix info_field1 info_field2 info_field3`
}

func (m *manipInfo) SetFlags(f *flag.FlagSet) {
	f.StringVar(&m.operator, "operator", "", "how to combine fields (max, min, mean)")
	f.StringVar(&m.prefix, "prefix", "", "prefix of new field being created")
}

func (m *manipInfo) Execute(_ context.Context, f *flag.FlagSet, _ ...interface{}) subcommands.ExitStatus {
	rdr, err := vcfgo.NewReader(os.Stdin, false)
	if err != nil {
		panic(err)
	}

	switch m.operator {
	case "max":
		rdr.AddInfoToHeader(m.prefix+"_max", "1", "Float", m.prefix+" max")
		rdr.AddInfoToHeader(m.prefix+"_max_name", "1", "String", "which "+m.prefix+" was the max")
	case "min":
		rdr.AddInfoToHeader(m.prefix+"_min", "1", "Float", m.prefix+" min")
		rdr.AddInfoToHeader(m.prefix+"_min_name", "1", "String", "which "+m.prefix+" was the min")
	case "mean":
		rdr.AddInfoToHeader(m.prefix+"_mean", "1", "Float", m.prefix+" mean")
		rdr.AddInfoToHeader(m.prefix+"_mean_name", "1", "String", "which "+m.prefix+" was the mean")
	default:
		return subcommands.ExitFailure
	}

	wrt, err := vcfgo.NewWriter(os.Stdout, rdr.Header)
	if err != nil {
		panic(err)
	}

	fields := f.Args()

	for {
		variant := rdr.Read()
		if variant == nil {
			break
		}

		var vals []float64
		var names []string

		for _, name := range fields {
			valI, _ := variant.Info().Get(name)
			val, ok := valI.(float64)
			if ok {
				vals = append(vals, val)
				names = append(names, name)
			}
		}

		var outVal float64
		var outName string

		if len(vals) != 0 {
			switch m.operator {
			case "max":
				outVal = vals[0]
				outName = names[0]
				for i, val := range vals {
					if val > outVal {
						outVal = val
						outName = names[i]
					}
				}
			case "min":
				outVal = vals[0]
				outName = names[0]
				for i, val := range vals {
					if val < outVal {
						outVal = val
						outName = names[i]
					}
				}
			case "mean":
				var total float64
				for _, val := range vals {
					total += val
				}
				outVal = total / float64(len(vals))
				outName = "mean"
			}

			variant.Info().Set(m.prefix+"_"+m.operator, outVal)
			variant.Info().Set(m.prefix+"_"+m.operator+"_name", outName)
			wrt.WriteVariant(variant)
		} else {
			wrt.WriteVariant(variant)
		}
	}
	return subcommands.ExitSuccess
}

type rank struct{}

func (*rank) Name() string { return "rank" }
func (*rank) Synopsis() string {
	return "create new info field with rank of variantMake new info field based off other fields"
}
func (*rank) Usage() string {
	return `rank riskGene1 riskGene2 riskGeneN`
}

func (r *rank) SetFlags(f *flag.FlagSet) {}

func isDmis(v *vcfgo.Variant) bool {
	revelI, _ := v.Info().Get("REVEL_score")
	revel, ok := revelI.(float64)
	if !ok {
		revel = 0.0
	}

	caddI, _ := v.Info().Get("CADD_phred")
	cadd, ok := caddI.(float64)
	if !ok {
		cadd = 0.0
	}

	del := cadd >= 25.0 || revel >= 0.5

	csqI, _ := v.Info().Get("vep_Consequence")
	csq, ok := csqI.(string)
	if !ok {
		csq = "."
	}

	if csq == "missense_variant" && del {
		return true
	}
	return false
}

func isLGD(v *vcfgo.Variant) bool {
	impactI, _ := v.Info().Get("vep_IMPACT")
	impact, ok := impactI.(string)
	if !ok {
		impact = "."
	}
	if impact == "HIGH" {
		return true
	}
	return false
}

func isConstrained(v *vcfgo.Variant) bool {
	pliI, _ := v.Info().Get("gnomAD_pLI")
	pli, ok := pliI.(float64)
	if !ok {
		pli = 0.0
	}

	if pli >= 0.5 {
		return true
	}
	return false
}

func getGnomAD(v *vcfgo.Variant) float64 {
	gnomadAFI, _ := v.Info().Get("eAF_popmax")
	gnomadAF, ok := gnomadAFI.(float64)
	if !ok {
		gnomadGAFI, _ := v.Info().Get("gAF_popmax")
		gnomadGAF, ok := gnomadGAFI.(float64)
		if !ok {
			return 0.0
		}
		return gnomadGAF
	}
	return gnomadAF
}

func isRare(v *vcfgo.Variant) bool {
	gnomadAF := getGnomAD(v)

	topMedI, _ := v.Info().Get("TOPMed_AF")
	topMed, ok := topMedI.(float64)
	if !ok {
		topMed = 0.0
	}

	if gnomadAF <= 0.0001 && topMed < 0.001 {
		return true
	}

	return false
}

func isSpliceDamage(v *vcfgo.Variant) bool {
	//impactI, _ := v.Info().Get("vep_IMPACT")
	//impact, ok := impactI.(string)
	//if !ok {
	//	impact = "."
	//}

	//if impact == "MODIFIER" || impact == "LOW" {
	maxI, _ := v.Info().Get("spliceAI_max")
	max, ok := maxI.(float64)
	if !ok {
		max = 0.0
	}

	if max >= 0.2 {
		return true
	}
	//}
	return false
}

func groupOne(v *vcfgo.Variant, riskGenes []string) bool {
	if isDmis(v) || isLGD(v) {
		if isRare(v) {
			geneI, _ := v.Info().Get("vep_SYMBOL")
			gene, ok := geneI.(string)
			if !ok {
				gene = "."
			}

			for _, rGene := range riskGenes {
				if gene == rGene {
					return true
				}
			}
		}
	}
	return false
}

func groupTwo(v *vcfgo.Variant) bool {
	if isLGD(v) {
		if isRare(v) {
			if isConstrained(v) {
				return true
			}
		}
	}
	return false
}

func groupTwoPointFive(v *vcfgo.Variant) bool {
	if isDmis(v) || isSpliceDamage(v) {
		if isRare(v) {
			if isConstrained(v) {
				return true
			}
		}
	}
	return false
}

func groupThree(v *vcfgo.Variant) bool {
	gnomadAF := getGnomAD(v)

	topMedI, _ := v.Info().Get("TOPMed_AF")
	topMed, ok := topMedI.(float64)
	if !ok {
		topMed = 0.0
	}

	if gnomadAF > 0.01 || topMed > 0.01 {
		return false
	}

	r := true
	recI, _ := v.Info().Get("recessive")
	if recI == nil {
		r = false
	}

	xr := true
	xrecI, _ := v.Info().Get("x_recessive")
	if xrecI == nil {
		xr = false
	}

	if r || xr {
		phomI, _ := v.Info().Get("phom")
		phom, ok := phomI.(float64)
		if !ok {
			phom = 1.0
		}
		if phom < 0.002 {
			return true
		}
	}

	return false
}

func groupFour(v *vcfgo.Variant) bool {
	if isDmis(v) || isLGD(v) || isSpliceDamage(v) {
		if isRare(v) {
			return true
		}
	}
	return false
}

func groupFive(v *vcfgo.Variant) bool {
	csqI, _ := v.Info().Get("vep_Consequence")
	csq, ok := csqI.(string)
	if !ok {
		csq = "."
	}

	if isLGD(v) || csq == "missense_variant" || isSpliceDamage(v) {
		gnomadAF := getGnomAD(v)

		if gnomadAF <= 0.001 && gnomadAF >= 0.0001 {
			return true
		}
	}

	return false
}

func groupFivePointFive(v *vcfgo.Variant) bool {
	dnvI, _ := v.Info().Get("denovo")
	_, ok := dnvI.(string)
	if ok {
		return true
	}
	hqdnvI, _ := v.Info().Get("hq_denovo")
	_, ok = hqdnvI.(string)
	if ok {
		return true
	}
	return false
}


func groupSix(v *vcfgo.Variant) bool {
	r := true
	recI, _ := v.Info().Get("recessive")
	if recI == nil {
		r = false
	}

	xr := true
	xrecI, _ := v.Info().Get("x_recessive")
	if xrecI == nil {
		xr = false
	}

	if r || xr {
		gnomadAF := getGnomAD(v)
		if gnomadAF < 0.01 {
			if isDmis(v) || isSpliceDamage(v) || isLGD(v) {
				return true
			}
		}

		phomI, _ := v.Info().Get("phom")
		phom, ok := phomI.(float64)
		if !ok {
			phom = 1.0
		}
		if phom < 0.05 && phom >= 0.002 {
			return true
		}
	}

	return false
}

func groupThreeCompHet(v *vcfgo.Variant) bool {
	chI, _ := v.Info().Get("slivar_comphet")
	if chI != nil {
		pchetI, _ := v.Info().Get("pchet")
		pchet, ok := pchetI.(float64)
		if !ok {
			pchet = 1.0
		}

		if pchet < 0.002 {
			return true
		}
	}

	return false
}

func groupSixCompHet(v *vcfgo.Variant) bool {
	chI, _ := v.Info().Get("slivar_comphet")
	if chI != nil {
		pchetI, _ := v.Info().Get("pchet")
		pchet, ok := pchetI.(float64)
		if !ok {
			pchet = 1.0
		}

		if pchet < 0.05 && pchet >= 0.002 {
			return true
		}


		gnomadAF := getGnomAD(v)
		if gnomadAF < 0.01 {
			if isDmis(v) || isSpliceDamage(v) || isLGD(v) {
				return true
			}
		}
	}

	return false
}

func (r *rank) Execute(_ context.Context, f *flag.FlagSet, _ ...interface{}) subcommands.ExitStatus {
	riskGenes := f.Args()

	rdr, err := vcfgo.NewReader(os.Stdin, false)
	if err != nil {
		fmt.Println(err)
		return subcommands.ExitFailure
	}

	rdr.AddInfoToHeader("rank", "1", "Float", "variant classifications")
	rdr.AddInfoToHeader("comphet_rank", "1", "Float", "variant classifications for half of compound het")

	wrt, err := vcfgo.NewWriter(os.Stdout, rdr.Header)
	if err != nil {
		fmt.Println(err)
		return subcommands.ExitFailure
	}

	for {
		variant := rdr.Read()
		if variant == nil {
			break
		}

		var rank float64
		switch {
		case groupOne(variant, riskGenes):
			rank = 1.0
		case groupTwo(variant):
			rank = 2.0
		case groupTwoPointFive(variant):
			rank = 2.5
		case groupThree(variant):
			rank = 3.0
		case groupFour(variant):
			rank = 4.0
		case groupFive(variant):
			rank = 5.0
		case groupFivePointFive(variant):
			rank = 5.5
		case groupSix(variant):
			rank = 6.0
		}

		if rank != 0.0 {
			variant.Info().Set("rank", rank)
		}

		var rankCompHet float64
		switch {
		case groupThreeCompHet(variant):
			rankCompHet = 3.0
		case groupSixCompHet(variant):
			rankCompHet = 6.0
		}

		if rankCompHet != 0.0 {
			variant.Info().Set("comphet_rank", rankCompHet)
		}

		wrt.WriteVariant(variant)
	}
	return subcommands.ExitSuccess
}

func getVarID(v *vcfgo.Variant) string {
	return v.Chromosome + "-" + strconv.Itoa(int(v.Pos)) + "-" + v.Reference + "-" + v.Alternate[0]
}

type filterCompHet struct {}

func (*filterCompHet) Name() string { return "filterCompHet" }
func (*filterCompHet) Synopsis() string {
	return "remove compound het pairs that dont both have ranks"
}
func (*filterCompHet) Usage() string {
	return `filterCompHet`
}

func (fch *filterCompHet) SetFlags(f *flag.FlagSet) {}

func (fch *filterCompHet) Execute(_ context.Context, f *flag.FlagSet, _ ...interface{}) subcommands.ExitStatus {

	rdr, err := vcfgo.NewReader(os.Stdin, false)
	if err != nil {
		fmt.Println(err)
		return subcommands.ExitFailure
	}

	wrt, err := vcfgo.NewWriter(os.Stdout, rdr.Header)
	if err != nil {
		fmt.Println(err)
		return subcommands.ExitFailure
	}

	type compHetVariant struct {
		variant *vcfgo.Variant
		slivarStr string
	}

	type compHetPair struct {
		ch1 compHetVariant
		ch2 compHetVariant
		paired bool
	}

	compHetMap := map[string]compHetPair{}
	for {
		variant := rdr.Read()
		if variant == nil {
			break
		}

		_, err := variant.Info().Get("comphet_rank")
		if err == nil {
			chI, err := variant.Info().Get("slivar_comphet")
			if err != nil {
				panic("should be a slivar compound het vcf, i.e. all variants should have info field 'slivar_comphet'")
			}
			if chString, ok := chI.(string); ok {
				chId := strings.Split(chString, "/")[2]
				if pair, ok := compHetMap[chId]; ok {
					pair.ch2 = compHetVariant{
						variant:   variant,
						slivarStr: chString,
					}
					pair.paired = true
					compHetMap[chId] = pair
				} else {
					compHetMap[chId] = compHetPair{
						ch1: compHetVariant{
							variant: variant,
							slivarStr: chString,
						},
						paired:   false,
					}
				}
			} else {
				if chSlice, ok := chI.([]string); ok {
					for _, chString := range chSlice {
						chId := strings.Split(chString, "/")[2]
						if pair, ok := compHetMap[chId]; ok {
							pair.ch2 = compHetVariant{
								variant:   variant,
								slivarStr: chString,
							}
							pair.paired = true
							compHetMap[chId] = pair
						} else {
							compHetMap[chId] = compHetPair{
								ch1: compHetVariant{
									variant: variant,
									slivarStr: chString,
								},
								paired:   false,
							}
						}
					}
				}
			}
		}
	}

	toWrite := map[string]*vcfgo.Variant{}
	for _, p := range compHetMap {
		if p.paired {
			v1Id := getVarID(p.ch1.variant)
			v2Id := getVarID(p.ch2.variant)
			if _, ok := toWrite[v1Id]; ok {
				oCHS, _ := toWrite[v1Id].Info().Get("slivar_comphet")
				toWrite[v1Id].Info().Set("slivar_comphet", oCHS.(string) + "," + p.ch1.slivarStr)
			} else {
				p.ch1.variant.Info().Set("slivar_comphet", p.ch1.slivarStr)
				toWrite[v1Id] = p.ch1.variant
			}
			if _, ok := toWrite[v2Id]; ok {
				oCHS, _ := toWrite[v2Id].Info().Get("slivar_comphet")
				toWrite[v2Id].Info().Set("slivar_comphet", oCHS.(string) + "," + p.ch2.slivarStr)
			} else {
				p.ch2.variant.Info().Set("slivar_comphet", p.ch2.slivarStr)
				toWrite[v2Id] = p.ch2.variant
			}
		}
	}

	for _, v := range toWrite {
		wrt.WriteVariant(v)
	}

	return subcommands.ExitSuccess
}

type popscores struct {
	chet *float64
	dom  *float64
	rec  *float64
}

type psap2vcf struct {
	txt     string
	proband string
	vcf string
}

func (*psap2vcf) Name() string { return "psap2vcf" }
func (*psap2vcf) Synopsis() string {
	return "convert psap report txt to vcf, with popscores in info field"
}
func (*psap2vcf) Usage() string {
	return `psap2vcf`
}

func (p *psap2vcf) SetFlags(f *flag.FlagSet) {
	f.StringVar(&p.txt, "txt", "", "txt file to extract psap values from")
	f.StringVar(&p.proband, "proband", "", "proband name")
	f.StringVar(&p.vcf, "vcf", "", "output vcf")
}

func (p *psap2vcf) Execute(_ context.Context, f *flag.FlagSet, _ ...interface{}) subcommands.ExitStatus {
	file, err := os.Open(p.txt)
	if err != nil {
		panic(err)
	}

	defer file.Close()

	psapM := make(map[string]*popscores)

	scanner := bufio.NewScanner(file)

	// find proband column
	var modelIdx int
	var scoreIdx int
	scanner.Scan()
	for idx, name := range strings.Split(scanner.Text(), "\t") {
		switch name {
		case "Dz.Model." + p.proband:
			modelIdx = idx
		case "popScore." + p.proband:
			scoreIdx = idx
		}
	}

	if modelIdx == 0 {
		panic("could not find proband column using supplied proband name")
	}

	for scanner.Scan() {
		ls := strings.Split(scanner.Text(), "\t")
		vID := ls[0] + "-" + ls[1] + "-" + ls[3] + "-" + ls[4]

		pType := ls[modelIdx]
		score := ls[scoreIdx]
		if _, ok := psapM[vID]; ok {
			switch pType {
			case "DOM-het":
				dom, err := strconv.ParseFloat(score, 64)
				if err != nil {
					panic(err)
				}
				psapM[vID].dom = &dom
			case "REC-hom":
				rec, err := strconv.ParseFloat(score, 64)
				if err != nil {
					panic(err)
				}
				psapM[vID].rec = &rec
			case "REC-chet":
				chet, err := strconv.ParseFloat(score, 64)
				if err != nil {
					panic(err)
				}
				psapM[vID].chet = &chet
			}
		} else {
			pops := popscores{}
			switch pType {
			case "DOM-het":
				dom, err := strconv.ParseFloat(score, 64)
				if err != nil {
					panic(err)
				}
				pops.dom = &dom
			case "REC-hom":
				rec, err := strconv.ParseFloat(score, 64)
				if err != nil {
					panic(err)
				}
				pops.rec = &rec
			case "REC-chet":
				chet, err := strconv.ParseFloat(score, 64)
				if err != nil {
					panic(err)
				}
				pops.chet = &chet
			}
			psapM[vID] = &pops
		}
	}

	hdr := vcfgo.NewHeader()
	hdr.FileFormat = "4.2"

	hdr.Infos["pdom"] = &vcfgo.Info{
		Id:          "pdom",
		Description: "psap dominate score",
		Number:      "1",
		Type:        "Float",
	}
	hdr.Infos["phom"] = &vcfgo.Info{
		Id:          "phom",
		Description: "psap homo score",
		Number:      "1",
		Type:        "Float",
	}
	hdr.Infos["pchet"] = &vcfgo.Info{
		Id:          "pchet",
		Description: "psap compound het score",
		Number:      "1",
		Type:        "Float",
	}

	wrt, err := vcfgo.NewWriter(os.Stdout, hdr)
	if err != nil {
		fmt.Println(err)
		return subcommands.ExitFailure
	}


	for k, v := range psapM {
		vID := strings.Split(k, "-")
		chrom := vID[0]
		pos, err := strconv.Atoi(vID[1])
		if err != nil {
			panic(err)
		}
		ref := vID[2]
		alt := vID[3]

		var variant *vcfgo.Variant

		variant = &vcfgo.Variant{
			Chromosome: chrom,
			Pos:        uint64(pos),
			Id_:        ".",
			Reference:  ref,
			Alternate:  []string{alt},
			Header: hdr,
			Filter: ".",
			Info_: vcfgo.NewInfoByte([]byte{}, hdr),
		}

		if v.dom != nil {
			_ = variant.Info().Set("pdom", *v.dom)
		}
		if v.rec != nil {
			_ = variant.Info().Set("phom", *v.rec)
		}
		if v.chet != nil {
			_ = variant.Info().Set("pchet", *v.chet)
		}
		wrt.WriteVariant(variant)
	}
	return subcommands.ExitSuccess
}

type coords struct {
	label string
}

func (*coords) Name() string { return "coords" }
func (*coords) Synopsis() string {
	return "add coordinates of variant to info field"
}
func (*coords) Usage() string {
	return `coords`
}

func (c *coords) SetFlags(f *flag.FlagSet) {
	f.StringVar(&c.label, "label", "", "label of coords i.e. hg19 -> hg19_pos")
}

func (c *coords) Execute(_ context.Context, f *flag.FlagSet, _ ...interface{}) subcommands.ExitStatus {
	rdr, err := vcfgo.NewReader(os.Stdin, false)
	if err != nil {
		fmt.Println(err)
		return subcommands.ExitFailure
	}

	rdr.AddInfoToHeader(c.label+"_chr", "1", "String", "chromosome from "+c.label)
	rdr.AddInfoToHeader(c.label+"_pos", "1", "Integer", "position from "+c.label)

	wrt, err := vcfgo.NewWriter(os.Stdout, rdr.Header)
	if err != nil {
		return subcommands.ExitFailure
	}

	for {
		variant := rdr.Read()
		if variant == nil {
			break
		}

		_ = variant.Info().Set(c.label+"_chr", variant.Chromosome)
		_ = variant.Info().Set(c.label+"_pos", int(variant.Pos))
		wrt.WriteVariant(variant)
	}
	return subcommands.ExitSuccess
}

type anchor struct {
	character string
	reference string
}

func (*anchor) Name() string { return "anchor" }
func (*anchor) Synopsis() string {
	return "remove charcter in vcf (*, -) and replace with anchored ref/alt"
}
func (*anchor) Usage() string {
	return `anchor -character "*" -reference ref.fa`
}

func (a *anchor) SetFlags(f *flag.FlagSet) {
	f.StringVar(&a.character, "character", "*", "character to replace")
	f.StringVar(&a.reference, "reference", "", "reference to get anchor base from")
}

func (a *anchor) Execute(_ context.Context, f *flag.FlagSet, _ ...interface{}) subcommands.ExitStatus {
	fa, err := faidx.New(a.reference)
	if err != nil {
		fmt.Println(err)
		return subcommands.ExitFailure
	}

	rdr, err := vcfgo.NewReader(os.Stdin, false)
	if err != nil {
		fmt.Println(err)
		return subcommands.ExitFailure
	}

	wrt, err := vcfgo.NewWriter(os.Stdout, rdr.Header)
	if err != nil {
		fmt.Println(err)
		return subcommands.ExitFailure
	}

	for {
		variant := rdr.Read()
		if variant == nil {
			break
		}
		switch a.character {
		case variant.Alt()[0]:
			bp, err := fa.Get(variant.Chromosome, int(variant.Pos)-2, int(variant.Pos)-1)
			if err != nil {
				fmt.Println(err)
				return subcommands.ExitFailure
			}

			variant.Pos = variant.Pos - 1
			variant.Reference = bp + variant.Ref()
			variant.Alternate = []string{bp}
			variant.Id_ = "."

			wrt.WriteVariant(variant)

		case variant.Reference:
			bp, err := fa.Get(variant.Chromosome, int(variant.Pos)-2, int(variant.Pos)-1)
			if err != nil {
				panic(err)
			}

			variant.Pos = variant.Pos - 1
			variant.Reference = bp
			variant.Alternate = []string{bp + variant.Alt()[0]}
			variant.Id_ = "."

			wrt.WriteVariant(variant)
		default:
			wrt.WriteVariant(variant)
		}
	}
	return subcommands.ExitSuccess
}

type mkVcf struct {
	//variants string
	//pedigree string
}

func (*mkVcf) Name() string { return "mkVcf" }
func (*mkVcf) Synopsis() string {
	return "take variants in the format 1-3453452-G-A-sampleId with optional pedigree file and outputs vcf"
}
func (*mkVcf) Usage() string {
	return `mkVcf -variants /path/to/variants.txt -pedigree /path/to/pedigree.ped`
}

func (v *mkVcf) SetFlags(f *flag.FlagSet) {
	//f.StringVar(&v.variants, "variants", "", "list of variants to convert")
	//f.StringVar(&v.pedigree, "pedigree", "", "pedigree file")
}

func (v *mkVcf) Execute(_ context.Context, f *flag.FlagSet, _ ...interface{}) subcommands.ExitStatus {
	hdr := vcfgo.NewHeader()
	hdr.FileFormat = "4.2"

	hdr.Infos["sample"] = &vcfgo.Info{
		Id:          "sample",
		Description: "samples",
		Number:      ".",
		Type:        "String",
	}

	wrt, err := vcfgo.NewWriter(os.Stdout, hdr)
	if err != nil {
		panic(err)
	}

	//type sv struct {
	//	chromosome string
	//	position int
	//	reference string
	//	alternate string
	//	sample string
	//}

	//var vl []sv

	scanner := bufio.NewScanner(os.Stdin)
	for scanner.Scan() {
		if strings.Contains(scanner.Text(), "#") {
			continue
		}

		ll := strings.Split(scanner.Text(), "-")
		chrom := ll[0]
		pos, err := strconv.Atoi(ll[1])
		if err != nil {
			log.Fatal("could not convert position string to int, something is wrong with format")
		}
		ref := ll[2]
		alt := ll[3]
		sample := ll[4]

		//cv := sv{
		//	chromosome: chrom,
		//	position:   pos,
		//	reference:  ref,
		//	alternate:  alt,
		//	sample:     sample,
		//}

		//vl = append(vl, cv)

		var variant *vcfgo.Variant
		variant = &vcfgo.Variant{
			Chromosome: chrom,
			Pos:        uint64(pos),
			Id_:        ".",
			Reference:  ref,
			Alternate:  []string{alt},
			Header:     hdr,
			Filter:     ".",
			Info_:      vcfgo.NewInfoByte([]byte{}, hdr),
		}
		_ = variant.Info().Set("sample", sample)
		wrt.WriteVariant(variant)
	}
	return subcommands.ExitSuccess
}

func mkCSQ(keys, vals []string) map[string]string {
	csq := map[string]string{}
	for i, k := range keys {
		csq[k] = vals[i]
	}
	return csq
}

func getCanon(csq []map[string]string) []map[string]string {
	rcsq := make([]map[string]string, 0)
	for _, c := range csq {
		if c["CANONICAL"] == "YES" {
			rcsq = append(rcsq, c)
		}
	}

	return rcsq
}

func getAppris(csq []map[string]string) []map[string]string {
	m := map[string]int{
		"P1": 7,
		"P2": 6,
		"P3": 5,
		"P4": 4,
		"P5": 3,
		"ALT1": 2,
		"ALT2": 1,
	}

	max := 0
	rcsq := make([]map[string]string, 0)
	for _, c := range csq {
		query := c["APPRIS"]
		if m[query] > max {
			rcsq = []map[string]string{c}
			max = m[query]
		}
		if m[query] == max {
			rcsq = append(rcsq, c)
		}
	}
	return rcsq
}

func getTSL(csq []map[string]string) []map[string]string {
	m := map[string]int{
		"1": 6,
		"2": 5,
		"3": 4,
		"4": 3,
		"5": 2,
		"NA": 1,
	}

	max := 0
	rcsq := make([]map[string]string, 0)
	for _, c := range csq {
		query := c["APPRIS"]
		if m[query] > max {
			rcsq = []map[string]string{c}
			max = m[query]
		}
		if m[query] == max {
			rcsq = append(rcsq, c)
		}
	}
	return rcsq
}

func getBiotype(csq []map[string]string) []map[string]string {
	rcsq := make([]map[string]string, 0)
	for _, c := range csq {
		if c["BIOTYPE"] == "protein_coding" {
			rcsq = append(rcsq, c)
		}
	}

	return rcsq
}

func getSevere(csq []map[string]string) []map[string]string {
	m := map[string]int{
		"transcript_ablation":                36,
		"splice_acceptor_variant":            35,
		"splice_donor_variant":               34,
		"stop_gained":                        33,
		"frameshift_variant":                 32,
		"stop_lost":                          31,
		"start_lost":                         30,
		"transcript_amplification":           29,
		"inframe_insertion":                  28,
		"inframe_deletion":                   27,
		"missense_variant":                   26,
		"protein_altering_variant":           25,
		"splice_region_variant":              24,
		"incomplete_terminal_codon_variant":  23,
		"start_retained_variant":             22,
		"stop_retained_variant":              21,
		"synonymous_variant":                 20,
		"coding_sequence_variant":            19,
		"mature_miRNA_variant":               18,
		"5_prime_UTR_variant":                17,
		"3_prime_UTR_variant":                16,
		"non_coding_transcript_exon_variant": 15,
		"intron_variant":                     14,
		"NMD_transcript_variant":             13,
		"non_coding_transcript_variant":      12,
		"upstream_gene_variant":              11,
		"downstream_gene_variant":            10,
		"TFBS_ablation":                      9,
		"TFBS_amplification":                 8,
		"TF_binding_site_variant":            7,
		"regulatory_region_ablation":         6,
		"regulatory_region_amplification":    5,
		"feature_elongation":                 4,
		"regulatory_region_variant":          3,
		"feature_truncation":                 2,
		"intergenic_variant":                 1,
	}

	max := 0
	rcsq := make([]map[string]string, 0)
	for _, c := range csq {
		query := strings.Split(c["Consequence"], "&")[0] //split multiple consequences, VEP always put most severe first
		if m[query] > max {
			rcsq = []map[string]string{c}
			max = m[query]
		}
		if m[query] == max {
			rcsq = append(rcsq, c)
		}
	}
	return rcsq
}

func rankCanon(csq []map[string]string) map[string]string {
	canons := getCanon(csq)
	if len(canons) == 1 {
		return canons[0]
	}
	if len(canons) == 0 {
		canons = csq
	}

	apps := getAppris(canons)
	if len(apps) == 1 {
		return apps[0]
	}

	tsl := getTSL(apps)
	if len(tsl) == 1 {
		return tsl[0]
	}

	biotype := getBiotype(tsl)
	if len(biotype) == 1 {
		return tsl[0]
	}
	if len(biotype) == 0 {
		biotype = tsl
	}

	severe := getSevere(biotype)
	return severe[0]
}

func rankSevere(csq []map[string]string) map[string]string {
	severe := getSevere(csq)
	if len(severe) == 1 {
		return severe[0]
	}

	canons := getCanon(severe)
	if len(canons) == 1 {
		return canons[0]
	}
	if len(canons) == 0 {
		canons = severe
	}

	apps := getAppris(canons)
	if len(apps) == 1 {
		return apps[0]
	}

	tsl := getTSL(apps)
	if len(tsl) == 1 {
		return tsl[0]
	}

	biotype := getBiotype(tsl)
	if len(biotype) >= 1 {
		return biotype[0]
	}
	return tsl[0]
}

type pullCSQ struct {
	extract string
}

func (*pullCSQ) Name() string { return "pullCSQ" }
func (*pullCSQ) Synopsis() string {
	return ""
}
func (*pullCSQ) Usage() string {
	return `pullCSQ -extract csqField1,csqField2`
}

func (p *pullCSQ) SetFlags(f *flag.FlagSet) {
	f.StringVar(&p.extract, "extract", "", "comma sep csq fields to extract")
}

func (p *pullCSQ) Execute(_ context.Context, f *flag.FlagSet, _ ...interface{}) subcommands.ExitStatus {
	// create VCF read, only read from stdin currently
	rdr, err := vcfgo.NewReader(os.Stdin, false)
	if err != nil {
		panic(err)
	}

	// parse fields from argument into array of fields to extract from csq
	extractFields := strings.Split(p.extract, ",")
	for _, f := range extractFields {
		rdr.AddInfoToHeader(f, "1", "String", "extracted from CSQ")
	}

	for _, f := range extractFields {
		rdr.AddInfoToHeader("canonical_" + f, "1", "String", "canonical " + f + " pulled from csq")
		rdr.AddInfoToHeader(f, "1", "String", "most severe " + f + " pulled from csq")
	}

	// create writer
	wrt, err := vcfgo.NewWriter(os.Stdout, rdr.Header)
	if err != nil {
		panic(err)
	}

	// get the csq key from the vcf header
	var csqInfo string
	if csqH, ok := rdr.Header.Infos["CSQ"]; ok {
		csqInfo = csqH.Description
	} else {
		log.Fatal("no CSQ field, please annotate with VEP")
	}
	csqKeys := strings.Split(strings.Split(csqInfo, "Format: ")[1], "|")


	for {
		variant := rdr.Read()
		if variant == nil {
			break
		}
		csq, err := variant.Info().Get("CSQ")
		if err != nil {
			wrt.WriteVariant(variant)
			continue
		}

		var scsq map[string]string
		var ccsq map[string]string

		switch csq := csq.(type) {
		case []string:
			acsq :=make([]map[string]string, 0)
			for _, c := range csq {
				acsq = append(acsq, mkCSQ(csqKeys, strings.Split(c, "|")))
			}
			scsq = rankSevere(acsq)
			ccsq = rankCanon(acsq)
		case string:
			scsq = mkCSQ(csqKeys, strings.Split(csq, "|"))
			ccsq = mkCSQ(csqKeys, strings.Split(csq, "|"))
		}

		for _, f := range extractFields {
			if ccsq[f] != "" {
				_ = variant.Info().Set("canonical_"+f, ccsq[f])
			}
			if scsq[f] != "" {
				_ = variant.Info().Set(f, scsq[f])
			}
		}
		wrt.WriteVariant(variant)
	}
	return subcommands.ExitSuccess
}

func main() {
	subcommands.Register(subcommands.HelpCommand(), "")
	subcommands.Register(subcommands.FlagsCommand(), "")
	subcommands.Register(subcommands.CommandsCommand(), "")
	subcommands.Register(&manipInfo{}, "")
	subcommands.Register(&rank{}, "")
	subcommands.Register(&anchor{}, "")
	subcommands.Register(&psap2vcf{}, "")
	subcommands.Register(&coords{}, "")
	subcommands.Register(&filterCompHet{}, "")
	subcommands.Register(&mkVcf{}, "")
	subcommands.Register(&pullCSQ{}, "")

	flag.Parse()
	ctx := context.Background()
	os.Exit(int(subcommands.Execute(ctx)))
}
