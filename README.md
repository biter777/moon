# moon

- Moon phase calculation
- Moon age calculation
- Angular diameter
- Distance to Earth
- Distance to Sun
- True Moon longitude
- Luminous Fraction of the Moon
- Physical data and parameters of the Moon
- Emoji supported
- Supported languages: English, French, German, Spanish, Italian, Russian, Ukrainian


[![GoDev](https://img.shields.io/badge/godev-reference-5b77b3)](https://pkg.go.dev/github.com/biter777/moon?tab=doc)
[![DOI](https://zenodo.org/badge/182808313.svg)](https://zenodo.org/badge/latestdoi/182808313)
[![GolangCI](https://lift.sonatype.com/api/badge/github.com/biter777/moon)](https://lift.sonatype.com/results/github.com/biter777/moon)
[![GolangCI](https://golangci.com/badges/github.com/biter777/moon.svg?style=flat)](https://golangci.com/r/github.com/biter777/moon)
[![GoReport](https://goreportcard.com/badge/github.com/biter777/moon)](https://goreportcard.com/report/github.com/biter777/moon)
[![Codiga](https://img.shields.io/badge/codiga%20quality-A+-brightgreen)](https://app.codiga.io/project/3255/dashboard)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/08eb1d2ff62e465091b3a288ae078a96)](https://www.codacy.com/manual/biter777/moon?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=biter777/moon&amp;utm_campaign=Badge_Grade)
[![License](https://img.shields.io/badge/License-BSD%202--Clause-brightgreen.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![Build status](https://ci.appveyor.com/api/projects/status/t9lpor9o8tpacpmr/branch/master?svg=true)](https://ci.appveyor.com/project/biter777/moon/branch/master)
[![CLDR](https://img.shields.io/badge/deepsource-passing-brightgreen)]([https://cldr.unicode.org/](https://deepsource.io/gh/biter777/moon))
<a href="//www.dmca.com/Protection/Status.aspx?ID=7a019cc5-ec73-464b-9707-4b33726f348f" title="DMCA.com Protection Status" class="dmca-badge"> <img src ="https://img.shields.io/badge/DMCA-protected-brightgreen" alt="DMCA.com Protection Status" /></a>
[![Dependencies Free](https://img.shields.io/badge/dependencies-free-brightgreen)](https://pkg.go.dev/github.com/biter777/moon?tab=imports)
[![Gluten Free](https://img.shields.io/badge/gluten-free-brightgreen)](https://www.scsglobalservices.com/services/gluten-free-certification)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen)](https://github.com/biter777/moon/pulls)
[![DepShield Badge](https://depshield.sonatype.org/badges/biter777/moon/depshield.svg)](https://depshield.github.io)
[![Stars](https://img.shields.io/github/stars/biter777/moon?label=Please%20like%20us&style=social)](https://github.com/biter777/moon/stargazers)
<br/>

## installation

```shell
go get github.com/biter777/moon
```
```go
import "github.com/biter777/moon"
```


## usage

```go
m := moon.New(time.Now())
fmt.Printf("Moon phase with emoji: %v\n", m.PhaseName(language.English, true))
fmt.Printf("Moon illumination: %v\n", int(m.Illumination*100)+"%")
fmt.Printf("Moon day (Moon age): %v\n", int(m.Age*100))
fmt.Printf("Full Moon in the current cycle at: %v\n", m.FullMoon())
fmt.Printf("Full Moon in the next cycle at: %v\n", m.NextFullMoon())
fmt.Printf("New Moon in the current cycle at: %v\n", m.NewMoon())
fmt.Printf("New Moon in the next cycle at: %v\n", m.NextNewMoon())
fmt.Printf("Diameter: %v\n", m.Diameter)
fmt.Printf("Distance: %v\n", m.Distance)
fmt.Printf("Longitude: %v\n", m.Longitude)
fmt.Printf("Sun distance: %v\n", m.SunDistance)
fmt.Printf("Sun angular diameter: %v\n", m.SunAngularDiameter)
fmt.Printf("Density: %v\n", moon.Density)
fmt.Printf("AverageDiameter: %v\n", moon.AverageDiameter)
fmt.Printf("AverageDistance: %v\n", moon.AverageDistance)
fmt.Printf("Gravity: %v\n", moon.Gravity)
fmt.Printf("Mass: %v\n", moon.Mass)
fmt.Printf("Average radius: %v\n", moon.AverageRadius)
fmt.Printf("Square: %v\n", moon.Square)
fmt.Printf("Synodic month: %v\n", moon.SynodicMonth)

```

For more complex options, consult the [documentation](https://pkg.go.dev/github.com/biter777/moon?tab=doc).

## Contributing

1. **Welcome pull requests, bug fixes and issue reports**

	[Contributors list](https://github.com/biter777/moon/graphs/contributors)
	
2. **Donate** - a donation isn't necessary, but it's welcome.

	<noscript><a href="https://liberapay.com/biter777/donate"><img alt="Donate using Liberapay" src="https://liberapay.com/assets/widgets/donate.svg"></a></noscript>
	[![ko-fi](https://www.ko-fi.com/img/githubbutton_sm.svg)](https://ko-fi.com/I2I61D1XZ) <a href="https://pay.cloudtips.ru/p/94fc4268" target="_blank"><img height="30" src="https://usa.visa.com/dam/VCOM/regional/lac/ENG/Default/Partner%20With%20Us/Payment%20Technology/visapos/full-color-800x450.jpg"></a> <a href="https://pay.cloudtips.ru/p/94fc4268" target="_blank"><img height="30" src="https://brand.mastercard.com/content/dam/mccom/brandcenter/thumbnails/mastercard_debit_sym_decal_web_105px.png"></a> <a href="https://pay.cloudtips.ru/p/94fc4268" target="_blank"><img height="30" src="https://developer.apple.com/assets/elements/icons/apple-pay/apple-pay.svg"></a> <a href="https://pay.cloudtips.ru/p/94fc4268" target="_blank"><img height="30" src="https://developers.google.com/pay/api/images/brand-guidelines/google-pay-mark.png"></a> <br/>

3. **Star us** - give us a star, please, if it's not against your religion :)


	[![Stars](https://img.shields.io/github/stars/biter777/moon?label=Please%20like%20us&style=social)](https://github.com/biter777/moon/stargazers)

