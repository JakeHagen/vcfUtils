package main

import (
	"context"
	"flag"
	"fmt"
	"os"

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
	AF, ok := gnomadAFI.(float64)
	if !ok {
		gnomadGAFI, _ := v.Info().Get("gAF_popmax")
		AF, ok = gnomadGAFI.(float64)
		if !ok {
			AF = 0.0
		}
	}
	return AF
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
	impactI, _ := v.Info().Get("vep_IMPACT")
	impact, ok := impactI.(string)
	if !ok {
		impact = "."
	}

	if impact == "MODIFIER" || impact == "LOW" {
		maxI, _ := v.Info().Get("spliceAI_max")
		max, ok := maxI.(float64)
		if !ok {
			max = 0.0
		}

		if max >= 0.5 {
			return true
		}
	}
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
	_, ok = recI.(string)
	if !ok {
		r = false
	}

	xr := true
	xrecI, _ := v.Info().Get("x_recessive")
	_, ok = xrecI.(string)
	if !ok {
		xr = false
	}

	ch := true
	chI, _ := v.Info().Get("slivar_comphet")
	_, ok = chI.(string)
	if !ok {
		ch = false
	}

	if r || ch || xr {
		if isLGD(v) {
			return true
		}

		pchetI, _ := v.Info().Get("pchet")
		pchet, ok := pchetI.(float64)
		if !ok {
			pchet = 1.0
		}

		if pchet <= 0.002 {
			return true
		}
		return false
	}
	return false
}

func groupFour(v *vcfgo.Variant) bool {
	if isDmis(v) || isLGD(v) {
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

	if isLGD(v) || csq == "missense" {
		gnomadAF := getGnomAD(v)

		if gnomadAF <= 0.001 && gnomadAF >= 0.0001 {
			return true
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

		var group float64
		switch {
		case groupOne(variant, riskGenes):
			group = 1.0
		case groupTwo(variant):
			group = 2.0
		case groupTwoPointFive(variant):
			group = 2.5
		case groupThree(variant):
			group = 3.0
		case groupFour(variant):
			group = 4.0
		case groupFive(variant):
			group = 5.0
		default:
			group = 6.0
		}

		variant.Info().Set("rank", group)

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

	flag.Parse()
	ctx := context.Background()
	os.Exit(int(subcommands.Execute(ctx)))
}
