package main

import (
	"bufio"
	"context"
	"flag"
	"fmt"
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
				variant.Info().Set("slivar_comphet", chString)
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

	flag.Parse()
	ctx := context.Background()
	os.Exit(int(subcommands.Execute(ctx)))
}
