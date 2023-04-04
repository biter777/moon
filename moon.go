// Package moon - package moon
package moon

import (
	"math"
	"time"

	"golang.org/x/text/language"
)

// Moon - Moon
type Moon struct {
	Age                float64 // age of the Moon, in days
	Diameter           float64 // angular diameter subtended by the Moon as seen by an observer at the centre of the Earth (degrees)
	Distance           float64 // distance of the Moon from the centre of the Earth (kilometres)
	Phase              float64 // Moon phase
	Illumination       float64 // illuminated fraction of the Moon (0 = New, 1 = Full)
	SunDistance        float64 // distance to Sun (kilometres)
	SunAngularDiameter float64 // angular diameter subtended by the Sun as seen by an observer at the centre of the Earth (degrees)
	Longitude          float64 // longitude
	pdata              float64
	timespace          float64
	quarters           [8]float64
}

// New - new instance
func New(t time.Time) *Moon {
	moon := &Moon{
		timespace: float64(t.Unix()),
		pdata:     utcToJulian(float64(t.Unix())),
	}

	// Calculation of the Sun's position
	day := moon.pdata - epoch                                       // Date within epoch
	n := fixAngle((360 / 365.2422) * day)                           // Mean anomaly of the Sun
	m := fixAngle(n + elonge - elongp)                              // Convert from perigee co-orginates to epoch 1980.0
	ec := kepler(m, eccent)                                         // Solve equation of Kepler
	ec = math.Sqrt((1+eccent)/(1-eccent)) * math.Tan(ec/2)          // ec
	ec = 2 * rad2deg(math.Atan(ec))                                 // True anomaly
	lambdasun := fixAngle(ec + elongp)                              // Sun's geocentric ecliptic longitude
	f := ((1 + eccent*math.Cos(deg2rad(ec))) / (1 - eccent*eccent)) // Orbital distance factor
	moon.SunDistance = sunsmax / f                                  // Distance to Sun in km
	moon.SunAngularDiameter = f * sunangsiz                         // Sun's angular size in degrees

	// Calsulation of the Moon's position
	ml := fixAngle(13.1763966*day + mmlong)               // Moon's mean longitude
	mm := fixAngle(ml - 0.1114041*day - mmlongp)          // Moon's mean anomaly
	ev := 1.2739 * math.Sin(deg2rad(2*(ml-lambdasun)-mm)) // Evection
	ae := 0.1858 * math.Sin(deg2rad(m))                   // Annual equation
	a3 := 0.37 * math.Sin(deg2rad(m))                     // Correction term
	mmP := mm + ev - ae - a3                              // Corrected anomaly
	mec := 6.2886 * math.Sin(deg2rad(mmP))                // Correction for the equation of the centre
	a4 := 0.214 * math.Sin(deg2rad(2*mmP))                // Another correction term
	lp := ml + ev + mec - ae + a4                         // Corrected longitude
	v := 0.6583 * math.Sin(deg2rad(2*(lp-lambdasun)))     // Variation
	moon.Longitude = lp + v                               // True longitude

	// Calculation of the phase of the Moon
	moonAge := moon.Longitude - lambdasun                                             // Age of the Moon in degrees
	moon.Illumination = (1 - math.Cos(deg2rad(moonAge))) / 2                          // Phase of the Moon, Illuminated fraction (0 to 1)
	moon.Distance = (msmax * (1 - mecc*mecc)) / (1 + mecc*math.Cos(deg2rad(mmP+mec))) // Distance of moon from the centre of the Earth, kilometres
	moonDFrac := moon.Distance / msmax                                                // moonDFrac
	moon.Diameter = mangsiz / moonDFrac                                               // Moon's angular diameter (degreees)
	moon.Phase = fixAngle(moonAge) / 360                                              // Phase (0 to 1)
	moon.Age = SynodicMonth * moon.Phase                                              // Age of moon (days)
	moon.phaseHunt()
	return moon
}

func rad2deg(r float64) float64 {
	return (r * 180) / math.Pi
}

func deg2rad(d float64) float64 {
	return (d * math.Pi) / 180
}

func fixAngle(a float64) float64 {
	return (a - 360*math.Floor(a/360))
}

func kepler(m, ecc float64) float64 {
	const epsilon = 0.000001
	m = deg2rad(m)
	e := m
	delta := e - ecc*math.Sin(e) - m
	e -= delta / (1 - ecc*math.Cos(e))
	for math.Abs(delta) > epsilon {
		delta = e - ecc*math.Sin(e) - m
		e -= delta / (1 - ecc*math.Cos(e))
	}
	return e
}

func (m *Moon) phaseHunt() {
	sdate := utcToJulian(m.timespace)
	adate := sdate - 45
	ats := m.timespace - 86400*45
	t := time.Unix(int64(ats), 0)
	yy := float64(t.Year())
	mm := float64(t.Month())

	k1 := math.Floor(float64(yy+((mm-1)*(1/12))-1900) * 12.3685)
	nt1 := meanPhase(adate, k1)
	adate = nt1

	var nt2, k2 float64
	for {
		adate += SynodicMonth
		k2 = k1 + 1
		nt2 = meanPhase(adate, k2)
		if math.Abs(nt2-sdate) < 0.75 {
			nt2 = truePhase(k2, 0.0)
		}
		if nt1 <= sdate && nt2 > sdate {
			break
		}
		nt1 = nt2
		k1 = k2
	}

	data := [8]float64{
		truePhase(k1, 0.0),
		truePhase(k1, 0.25),
		truePhase(k1, 0.5),
		truePhase(k1, 0.75),
		truePhase(k2, 0.0),
		truePhase(k2, 0.25),
		truePhase(k2, 0.5),
		truePhase(k2, 0.75),
	}

	for i := 0; i < 8; i++ {
		m.quarters[i] = (data[i] - 2440587.5) * 86400 // convert to UNIX time
	}
}

func utcToJulian(t float64) float64 {
	return t/86400 + 2440587.5
}

func julianToUtc(t float64) float64 {
	return t*86400 + 2440587.5
}

/*
Calculates time of the mean new Moon for a given
base date. This argument K to this function is the
precomputed synodic month index, given by:

	K = (year - 1900) * 12.3685

where year is expressed as a year aand fractional year
*/
func meanPhase(sdate float64, k float64) float64 {
	// Time in Julian centuries from 1900 January 0.5
	t := (sdate - 2415020.0) / 36525
	t2 := t * t
	t3 := t2 * t

	return 2415020.75933 + SynodicMonth*k +
		0.0001178*t2 -
		0.000000155*t3 +
		0.00033*math.Sin(deg2rad(166.56+132.87*t-0.009173*t2))
}

func truePhase(k float64, phase float64) float64 {
	k += phase       // Add phase to new moon time
	t := k / 1236.85 // Time in Julian centures from 1900 January 0.5
	t2 := t * t
	t3 := t2 * t
	pt := 2415020.75933 + SynodicMonth*k +
		0.0001178*t2 -
		0.000000155*t3 +
		0.00033*math.Sin(deg2rad(166.56+132.87*t-0.009173*t2))

	m := 359.2242 + 29.10535608*k - 0.0000333*t2 - 0.00000347*t3       // Sun's mean anomaly
	mprime := 306.0253 + 385.81691806*k + 0.0107306*t2 + 0.00001236*t3 // Moon's mean anomaly
	f := 21.2964 + 390.67050646*k - 0.0016528*t2 - 0.00000239*t3       // Moon's argument of latitude

	if phase < 0.01 || math.Abs(phase-0.5) < 0.01 {
		// Corrections for New and Full Moon
		pt += (0.1734-0.000393*t)*math.Sin(deg2rad(m)) +
			0.0021*math.Sin(deg2rad(2*m)) -
			0.4068*math.Sin(deg2rad(mprime)) +
			0.0161*math.Sin(deg2rad(2*mprime)) -
			0.0004*math.Sin(deg2rad(3*mprime)) +
			0.0104*math.Sin(deg2rad(2*f)) -
			0.0051*math.Sin(deg2rad(m+mprime)) -
			0.0074*math.Sin(deg2rad(m-mprime)) +
			0.0004*math.Sin(deg2rad(2*f+m)) -
			0.0004*math.Sin(deg2rad(2*f-m)) -
			0.0006*math.Sin(deg2rad(2*f+mprime)) +
			0.0010*math.Sin(deg2rad(2*f-mprime)) +
			0.0005*math.Sin(deg2rad(m+2*mprime))
	} else if math.Abs(phase-0.25) < 0.01 || math.Abs(phase-0.75) < 0.01 {
		pt += (0.1721-0.0004*t)*math.Sin(deg2rad(m)) +
			0.0021*math.Sin(deg2rad(2*m)) -
			0.6280*math.Sin(deg2rad(mprime)) +
			0.0089*math.Sin(deg2rad(2*mprime)) -
			0.0004*math.Sin(deg2rad(3*mprime)) +
			0.0079*math.Sin(deg2rad(2*f)) -
			0.0119*math.Sin(deg2rad(m+mprime)) -
			0.0047*math.Sin(deg2rad(m-mprime)) +
			0.0003*math.Sin(deg2rad(2*f+m)) -
			0.0004*math.Sin(deg2rad(2*f-m)) -
			0.0006*math.Sin(deg2rad(2*f+mprime)) +
			0.0021*math.Sin(deg2rad(2*f-mprime)) +
			0.0003*math.Sin(deg2rad(m+2*mprime)) +
			0.0004*math.Sin(deg2rad(m-2*mprime)) -
			0.0003*math.Sin(deg2rad(2*m+mprime))
		if phase < 0.5 { // First quarter correction
			pt += 0.0028 - 0.0004*math.Cos(deg2rad(m)) + 0.0003*math.Cos(deg2rad(mprime))
		} else { // Last quarter correction
			pt += -0.0028 + 0.0004*math.Cos(deg2rad(m)) - 0.0003*math.Cos(deg2rad(mprime))
		}
	}

	return pt
}

// NewMoon - the time of the New Moon in the current lunar cycle, i.e., the start of the current cycle
func (m *Moon) NewMoon() time.Time {
	return time.Unix(int64(m.quarters[0]), 0)
}

// FirstQuarter - the time of the first quarter in the current lunar cycle
func (m *Moon) FirstQuarter() time.Time {
	return time.Unix(int64(m.quarters[1]), 0)
}

// FullMoon - the time of the Full Moon in the current lunar cycle
func (m *Moon) FullMoon() time.Time {
	return time.Unix(int64(m.quarters[2]), 0)
}

// LastQuarter - the time of the last quarter in the current lunar cycle
func (m *Moon) LastQuarter() time.Time {
	return time.Unix(int64(m.quarters[3]), 0)
}

// NextNewMoon - the time of the New Moon in the next lunar cycle, i.e., the start of the next cycle
func (m *Moon) NextNewMoon() time.Time {
	return time.Unix(int64(m.quarters[4]), 0)
}

// NextFirstQuarter - the time of the first quarter in the next lunar cycle
func (m *Moon) NextFirstQuarter() time.Time {
	return time.Unix(int64(m.quarters[1]), 0)
}

// NextFullMoon - the time of the Full Moon in the next lunar cycle
func (m *Moon) NextFullMoon() time.Time {
	return time.Unix(int64(m.quarters[6]), 0)
}

// NextLastQuarter - the time of the last quarter in the next lunar cycle
func (m *Moon) NextLastQuarter() time.Time {
	return time.Unix(int64(m.quarters[7]), 0)
}

// PhaseName - the phase name
func (m *Moon) PhaseName(lang language.Tag, emoji bool) (name string) {
	i := int(math.Floor((m.Phase + 0.0625) * 8))
	switch lang {
	case language.Dutch:
		name = m.phaseNameDe(i)
	case language.Spanish, language.EuropeanSpanish:
		name = m.phaseNameEs(i)
	case language.French:
		name = m.phaseNameFr(i)
	case language.Italian:
		name = m.phaseNameIt(i)
	case language.Russian:
		name = m.phaseNameRu(i)
	case language.Ukrainian:
		name = m.phaseNameUa(i)
	default:
		name = m.phaseNameEn(i)
	}

	if emoji {
		emojs := [9]string{
			0: "🌑",
			1: "🌒",
			2: "🌓",
			3: "🌔",
			4: "🌕",
			5: "🌖",
			6: "🌗",
			7: "🌘",
			8: "🌑",
		}
		return emojs[i] + " " + name
	}

	return name
}

func (m *Moon) phaseNameEn(i int) string {
	names := [9]string{
		0: "New Moon",
		1: "Waxing Crescent",
		2: "First Quarter",
		3: "Waxing Gibbous",
		4: "Full Moon",
		5: "Waning Gibbous",
		6: "Third Quarter",
		7: "Waning Crescent",
		8: "New Moon",
	}
	return names[i]
}

func (m *Moon) phaseNameEs(i int) string {
	names := [9]string{
		0: "Luna Nueva",
		1: "Creciente Creciente",
		2: "Primer cuarto",
		3: "Cera menguante",
		4: "Luna llena",
		5: "Menguante menguante",
		6: "Tercer cuarto",
		7: "Creciente menguante",
		8: "Luna Nueva",
	}
	return names[i]
}

func (m *Moon) phaseNameDe(i int) string {
	names := [9]string{
		0: "Neumond",
		1: "Zunehmende Sichel",
		2: "Erstes Viertel",
		3: "Zunehmender Halbmond",
		4: "Vollmond",
		5: "Abnehmender Halbmond",
		6: "Letztes Viertel",
		7: "Abnehmende Sichel",
		8: "Neumond",
	}
	return names[i]
}

func (m *Moon) phaseNameFr(i int) string {
	names := [9]string{
		0: "Nouvelle Lune",
		1: "Croissant croissant",
		2: "Premier quart",
		3: "Gibbeuse croissante",
		4: "Pleine Lune",
		5: "Gibbeuse décroissante",
		6: "Troisième quart",
		7: "Croissant déclinant",
		8: "Nouvelle Lune",
	}
	return names[i]
}

func (m *Moon) phaseNameIt(i int) string {
	names := [9]string{
		0: "Luna Nuova",
		1: "Mezzaluna crescente",
		2: "Primo quarto",
		3: "Gibboso crescente",
		4: "Luna piena",
		5: "Gibboso calante",
		6: "Terzo quarto",
		7: "Mezzaluna calante",
		8: "Luna Nuova",
	}
	return names[i]
}

func (m *Moon) phaseNameRu(i int) string {
	names := [9]string{
		0: "Новолуние",
		1: "Растущая луна",
		2: "Первая четверть",
		3: "Растущая луна",
		4: "Полнолуние",
		5: "Убывающая луна",
		6: "Третья четверть",
		7: "Убывающая Луна",
		8: "Новолуние",
	}
	return names[i]
}

func (m *Moon) phaseNameUa(i int) string {
	names := [9]string{
		0: "Молодий місяць",
		1: "Зростаючий півмісяць",
		2: "Перша чверть",
		3: "Зростаючий гіббоус",
		4: "Повний місяць",
		5: "Waning Gibbous",
		6: "Третя чверть",
		7: "Спадаючий півмісяць",
		8: "Молодий місяць",
	}
	return names[i]
}

// ZodiacEmoji - zodiac emoji
func (m *Moon) ZodiacEmoji() string {
	return zodiacEmoji(m.getZodiacNumber())
}

func zodiacEmoji(i int) string {
	emojs := [13]string{
		0:  "♈", // Aries
		1:  "♉", // Taurus
		2:  "♊", // Gemini
		3:  "♋", // Cancer
		4:  "♌", // Leo
		5:  "♍", // Virgo
		6:  "♎", // Libra
		7:  "♏", // Scorpio
		8:  "♐", // Sagittarius
		9:  "♑", // Capricorn
		10: "♒", // Aquarius
		11: "♓", // Pisces
		12: "♈", // Aries
	}
	return emojs[i]
}

// Zodiac - zodiac sign
func (m *Moon) Zodiac(lang language.Tag, emoji bool) (zodiac string) {
	i := m.getZodiacNumber()
	switch lang {
	case language.French:
		zodiac = m.zodiacFr(i)
	case language.Dutch:
		zodiac = m.zodiacDe(i)
	case language.Spanish, language.EuropeanSpanish:
		zodiac = m.zodiacEs(i)
	case language.Italian:
		zodiac = m.zodiacIt(i)
	case language.Russian:
		zodiac = m.zodiacRu(i)
	case language.Ukrainian:
		zodiac = m.zodiacRu(i)
	default:
		zodiac = m.zodiacEn(i)
	}

	if emoji {
		return zodiacEmoji(i) + " " + zodiac
	}

	return zodiac
}

func (m *Moon) zodiacEn(i int) string {
	signs := [13]string{
		"Aries",
		"Taurus",
		"Gemini",
		"Cancer",
		"Leo",
		"Virgo",
		"Libra",
		"Scorpio",
		"Sagittarius",
		"Capricorn",
		"Aquarius",
		"Pisces",
		"Aries",
	}
	return signs[i]
}

func (m *Moon) zodiacFr(i int) string {
	signs := [13]string{
		"Bélier",
		"Taureau",
		"Gémeaux",
		"Cancer",
		"Leo",
		"Vierge",
		"Balance",
		"Scorpion",
		"Sagittaire",
		"Capricorne",
		"Verseau",
		"Poissons",
		"Bélier",
	}
	return signs[i]
}

func (m *Moon) zodiacDe(i int) string {
	signs := [13]string{
		"Widder",
		"Stier",
		"Zwillinge",
		"Krebs",
		"Löwe",
		"Jungfrau",
		"Waage",
		"Skorpion",
		"Schütze",
		"Steinbock",
		"Wassermann",
		"Fische",
		"Widder",
	}
	return signs[i]
}

func (m *Moon) zodiacEs(i int) string {
	signs := [13]string{
		"Aries",
		"Tauro",
		"Geminis",
		"Cáncer",
		"León",
		"Virgo",
		"Libra",
		"Escorpión",
		"Sagitario",
		"Capricornio",
		"Acuario",
		"Piscis",
		"Aries",
	}
	return signs[i]
}

func (m *Moon) zodiacIt(i int) string {
	signs := [13]string{
		"Ariete",
		"Toro",
		"Gemelli",
		"Cancro",
		"Leo",
		"Vergine",
		"Libra",
		"Scorpione",
		"Sagittario",
		"Capricorno",
		"Acquario",
		"Pesci",
		"Ariete",
	}
	return signs[i]
}

func (m *Moon) zodiacRu(i int) string {
	signs := [13]string{
		"Овен",
		"Телец",
		"Близнецы",
		"Рак",
		"Лев",
		"Дева",
		"Весы",
		"Скорпион",
		"Стрелец",
		"Козерог",
		"Водолей",
		"Рыбы",
		"Овен",
	}
	return signs[i]
}

func (m *Moon) zodiacUa(i int) string {
	signs := [13]string{
		"Овен",
		"Телець",
		"Близнюки",
		"Рак",
		"Лев",
		"Діва",
		"Терези",
		"Скорпіон",
		"Стрілець",
		"Козеріг",
		"Водолій",
		"Риби",
		"Овен",
	}
	return signs[i]
}

func (m *Moon) getZodiacNumber() int {
	if m.Longitude < 33.18 {
		return 0
	} else if m.Longitude < 51.16 {
		return 1
	} else if m.Longitude < 93.44 {
		return 2
	} else if m.Longitude < 119.48 {
		return 3
	} else if m.Longitude < 135.30 {
		return 4
	} else if m.Longitude < 173.34 {
		return 5
	} else if m.Longitude < 224.17 {
		return 6
	} else if m.Longitude < 242.57 {
		return 7
	} else if m.Longitude < 271.26 {
		return 8
	} else if m.Longitude < 302.49 {
		return 9
	} else if m.Longitude < 311.72 {
		return 10
	} else if m.Longitude < 348.58 {
		return 11
	} else {
		return 12
	}
}

// Positive - positive Moon moments
func (m *Moon) Positive(lang language.Tag) []string {
	switch lang {
	case language.Dutch:
		return m.positiveDe()
	case language.Spanish, language.EuropeanSpanish:
		return m.positiveEs()
	case language.French:
		return m.positiveFr()
	case language.Italian:
		return m.positiveIt()
	case language.Russian:
		return m.positiveRu()
	case language.Ukrainian:
		return m.positiveUa()
	default:
		return m.positiveEn()
	}
}

// Negative - negative Moon moments
func (m *Moon) Negative(lang language.Tag) []string {
	switch lang {
	case language.Dutch:
		return m.negativeDe()
	case language.Spanish, language.EuropeanSpanish:
		return m.negativeEs()
	case language.French:
		return m.negativeFr()
	case language.Italian:
		return m.negativeIt()
	case language.Russian:
		return m.negativeRu()
	case language.Ukrainian:
		return m.negativeUa()
	default:
		return m.negativeEn()
	}
}

func (m *Moon) positiveRu() []string {
	switch int(m.Age) {
	case 1:
		return []string{"путешествия", "зачатие"}
	case 2:
		return []string{"начинания", "финансы", "обучение", "творчество", "путешествия", "зачатие", "дом"}
	case 3:
		return []string{"нагрузки", "сад", "секс"}
	case 4:
		return []string{"дом", "финансы", "отдых", "зачатие", "сад"}
	case 5:
		return []string{"алко", "стрижка", "нагрузки", "сад"}
	case 6:
		return []string{"общение", "обучение", "творчество", "дом"}
	case 7:
		return []string{"общение", "сны", "дом", "отдых", "творчество"}
	case 8:
		return []string{"путешествия", "недвижимость", "отдых", "начинания", "творчество", "сны"}
	case 9:
		return []string{"нагрузки", "сад", "путешествия"}
	case 10:
		return []string{"начинания", "отдых", "дом", "финансы", "общение", "творчество", "зачатие", "недвижимость", "нагрузки", "сад"}
	case 11:
		return []string{"путешествия", "отдых", "стрижка", "сад", "зачатие", "сны", "творчество"}
	case 12:
		return []string{"сад", "сны", "стрижка", "общение"}
	case 13:
		return []string{"общение", "дом", "секс", "творчество", "обучение", "сны", "финансы", "алко", "стрижка"}
	case 14:
		return []string{"начинания", "финансы", "творчество", "путешествия", "дом", "общение", "обучение", "сад", "нагрузки", "зачатие"}
	case 15:
		return []string{"дом", "сад"}
	case 16:
		return []string{"общение", "отдых", "дом", "путешествия", "творчество", "брак", "зачатие", "сны"}
	case 17:
		return []string{"общение", "отдых", "секс", "алко", "творчество", "финансы", "недвижимость", "путешествия", "брак", "сны", "стрижка", "сад"}
	case 18:
		return []string{"путешествия", "творчество", "нагрузки", "сны"}
	case 19:
		return []string{"путешествия", "творчество", "дом", "сад", "сны"}
	case 20:
		return []string{"начинания", "общение", "финансы", "недвижимость", "творчество", "дом", "путешествия", "обучение", "сад"}
	case 21:
		return []string{"общение", "недвижимость", "путешествия", "отдых", "дом", "нагрузки", "начинания", "финансы", "секс", "брак", "зачатие", "сад"}
	case 22:
		return []string{"обучение", "алко", "зачатие", "нагрузки", "сны", "сад"}
	case 23:
		return []string{"дом", "сад"}
	case 24:
		return []string{"секс", "нагрузки", "общение", "финансы", "недвижимость", "творчество", "зачатие", "дом", "сад"}
	case 25:
		return []string{"недвижимость", "отдых", "зачатие", "сад"}
	case 26:
		return []string{"обучение", "брак", "зачатие", "дом", "сад", "сны"}
	case 27:
		return []string{"начинания", "финансы", "стрижка", "финансы", "путешествия", "брак", "сад"}
	case 28:
		return []string{"творчество", "зачатие", "сны", "начинания", "финансы", "недвижимость", "дом", "нагрузки"}
	case 29:
		return []string{"дом", "стрижка"}
	case 30:
		return []string{"отдых", "творчество", "обучение", "зачатие"}
	default:
		panic("moon age error")
	}
}

func (m *Moon) positiveUa() []string {
	switch int(m.Age) {
	case 1:
		return []string{"подорожі", "зачаття"}
	case 2:
		return []string{"починання", "фінанси", "навчання", "творчість", "подорожі", "зачаття", "будинок"}
	case 3:
		return []string{"навантаження", "сад", "секс"}
	case 4:
		return []string{"будинок", "фінанси", "відпочинок", "зачаття", "сад"}
	case 5:
		return []string{"алко", "стрижка", "навантаження", "сад"}
	case 6:
		return []string{"спілкування", "навчання", "творчість", "дім"}
	case 7:
		return []string{"спілкування", "сни", "будинок", "відпочинок", "творчість"}
	case 8:
		return []string{"подорожі", "нерухомість", "відпочинок", "починання", "творчість", "сни"}
	case 9:
		return []string{"навантаження", "сад", "подорожі"}
	case 10:
		return []string{"починання", "відпочинок", "будинок", "фінанси", "спілкування", "творчість", "зачаття", "нерухомість", "навантаження", "сад"}
	case 11:
		return []string{"подорожі", "відпочинок", "стрижка", "сад", "зачаття", "сни", "творчість"}
	case 12:
		return []string{"садок", "сни", "стрижка", "спілкування"}
	case 13:
		return []string{"спілкування", "дім", "секс", "творчість", "навчання", "сни", "фінанси", "алко", "стрижка"}
	case 14:
		return []string{"починання", "фінанси", "творчість", "подорожі", "будинок", "спілкування", "навчання", "сад", "навантаження", "зачаття"}
	case 15:
		return []string{"будинок", "сад"}
	case 16:
		return []string{"спілкування", "відпочинок", "будинок", "подорожі", "творчість", "шлюб", "зачаття", "сни"}
	case 17:
		return []string{"спілкування", "відпочинок", "секс", "алко", "творчість", "фінанси", "нерухомість", "подорожі", "шлюб", "сни", "стрижка", " сад"}
	case 18:
		return []string{"подорожі", "творчість", "навантаження", "сни"}
	case 19:
		return []string{"подорожі", "творчість", "будинок", "сад", "сни"}
	case 20:
		return []string{"починання", "спілкування", "фінанси", "нерухомість", "творчість", "будинок", "подорожі", "навчання", "сад"}
	case 21:
		return []string{"спілкування", "нерухомість", "подорожі", "відпочинок", "будинок", "навантаження", "починання", "фінанси", "секс", "шлюб", "зачаття", " сад"}
	case 22:
		return []string{"навчання", "алко", "зачаття", "навантаження", "сни", "сад"}
	case 23:
		return []string{"будинок", "сад"}
	case 24:
		return []string{"секс", "навантаження", "спілкування", "фінанси", "нерухомість", "творчість", "зачаття", "будинок", "сад"}
	case 25:
		return []string{"нерухомість", "відпочинок", "зачаття", "сад"}
	case 26:
		return []string{"навчання", "шлюб", "зачаття", "будинок", "сад", "сни"}
	case 27:
		return []string{"починання", "фінанси", "стрижка", "фінанси", "подорожі", "шлюб", "сад"}
	case 28:
		return []string{"творчість", "зачаття", "сни", "починання", "фінанси", "нерухомість", "будинок", "навантаження"}
	case 29:
		return []string{"будинок", "стрижка"}
	case 30:
		return []string{"відпочинок", "творчість", "навчання", "зачаття"}
	default:
		panic("moon age error")
	}
}

func (m *Moon) positiveDe() []string {
	switch int(m.Age) {
	case 1:
		return []string{"reise", "konzeption"}
	case 2:
		return []string{"Anfänge", "Finanzen", "Lernen", "Kreativität", "Reisen", "Konzeption", "Zuhause"}
	case 3:
		return []string{"load", "garden", "sex"}
	case 4:
		return []string{"Zuhause", "Finanzen", "Urlaub", "Empfängnis", "Garten"}
	case 5:
		return []string{"alco", "haircut", "load", "garden"}
	case 6:
		return []string{"Kommunikation", "Lernen", "Kreativität", "Zuhause"}
	case 7:
		return []string{"Kommunikation", "Träume", "Zuhause", "Ruhe", "Kreativität"}
	case 8:
		return []string{"travel", "Immobilien", "Urlaub", "Anfang", "Kreativität", "Träume"}
	case 9:
		return []string{"loads", "garden", "travel"}
	case 10:
		return []string{"Anfänge", "Urlaub", "Zuhause", "Finanzen", "Kommunikation", "Kreativität", "Konzeption", "Immobilien", "Belastung", "Garten"}
	case 11:
		return []string{"reise", "urlaub", "haarschnitt", "garten", "konzeption", "träume", "kreativität"}
	case 12:
		return []string{"garden", "dreams", "haircut", "socialization"}
	case 13:
		return []string{"Kommunikation", "Zuhause", "Sex", "Kreativität", "Lernen", "Träume", "Finanzen", "Alkohol", "Haarschnitt"}
	case 14:
		return []string{"Anfänge", "Finanzen", "Kreativität", "Reisen", "Zuhause", "Kommunikation", "Lernen", "Garten", "Belastung", "Konzeption"}
	case 15:
		return []string{"house", "garden"}
	case 16:
		return []string{"Kommunikation", "Urlaub", "Zuhause", "Reise", "Kreativität", "Ehe", "Empfängnis", "Träume"}
	case 17:
		return []string{"Kommunikation", "Urlaub", "Sex", "Alko", "Kreativität", "Finanzen", "Immobilien", "Reisen", "Ehe", "Träume", "Haarschnitt", "Garten"}
	case 18:
		return []string{"Reisen", "Kreativität", "Lasten", "Träume"}
	case 19:
		return []string{"travel", "art", "home", "garden", "dreams"}
	case 20:
		return []string{"Anfänge", "Kommunikation", "Finanzen", "Immobilien", "Kreativität", "Zuhause", "Reisen", "Lernen", "Garten"}
	case 21:
		return []string{"Kommunikation", "Immobilien", "Reise", "Urlaub", "Zuhause", "Arbeitsbelastung", "Startups", "Finanzen", "Sex", "Ehe", "Empfängnis", "Garten"}
	case 22:
		return []string{"training", "alco", "konzeption", "last", "träume", "garten"}
	case 23:
		return []string{"house", "garden"}
	case 24:
		return []string{"sex", "load", "communication", "finance", "real Estate", "creativity", "konzeption", "home", "garden"}
	case 25:
		return []string{"immobilie", "urlaub", "konzeption", "garten"}
	case 26:
		return []string{"lernen", "heiraten", "empfangen", "heimat", "garten", "träume"}
	case 27:
		return []string{"Anfänge", "Finanzen", "Haarschnitt", "Finanzen", "Reisen", "Ehe", "Garten"}
	case 28:
		return []string{"kreativität", "konzeption", "träume", "anfänge", "finanzen", "immobilien", "home", "lasten"}
	case 29:
		return []string{"house", "haircut"}
	case 30:
		return []string{"rest", "creativity", "learning", "konzeption"}
	default:
		panic("moon age error")
	}
}

func (m *Moon) positiveEs() []string {
	switch int(m.Age) {
	case 1:
		return []string{"viaje", "concepción"}
	case 2:
		return []string{"comienzos", "finanzas", "aprendizaje", "creatividad", "viajes", "concebir", "casa"}
	case 3:
		return []string{"cargar", "jardín", "sexo"}
	case 4:
		return []string{"casa", "finanzas", "vacaciones", "concepción", "jardín"}
	case 5:
		return []string{"alco", "corte de pelo", "carga", "jardín"}
	case 6:
		return []string{"comunicación", "aprendizaje", "creatividad", "inicio"}
	case 7:
		return []string{"comunicación", "sueños", "casa", "descanso", "creatividad"}
	case 8:
		return []string{"viaje", "bienes raíces", "vacaciones", "inicio", "creatividad", "sueños"}
	case 9:
		return []string{"cargas", "jardín", "viaje"}
	case 10:
		return []string{"comienzos", "vacaciones", "casa", "finanzas", "comunicación", "creatividad", "concepción", "bienes raíces", "carga", "jardín"}
	case 11:
		return []string{"viaje", "vacaciones", "corte de pelo", "jardín", "concepción", "sueños", "creatividad"}
	case 12:
		return []string{"jardín", "sueños", "corte de pelo", "socialización"}
	case 13:
		return []string{"comunicación", "hogar", "sexo", "creatividad", "aprendizaje", "sueños", "finanzas", "alcohol", "corte de pelo"}
	case 14:
		return []string{"comienzos", "finanzas", "creatividad", "viajes", "casa", "comunicación", "aprendizaje", "jardín", "carga", "concepción"}
	case 15:
		return []string{"casa", "jardín"}
	case 16:
		return []string{"comunicación", "vacaciones", "casa", "viaje", "creatividad", "matrimonio", "concepción", "sueños"}
	case 17:
		return []string{"comunicación", "vacaciones", "sexo", "alcohólico", "creatividad", "finanzas", "bienes raíces", "viajes", "matrimonio", "sueños", "corte de pelo", "jardín"}
	case 18:
		return []string{"viajes", "creatividad", "cargas", "sueños"}
	case 19:
		return []string{"viaje", "arte", "casa", "jardín", "sueños"}
	case 20:
		return []string{"comienzos", "comunicación", "finanzas", "bienes raíces", "creatividad", "hogar", "viajes", "aprendizaje", "jardín"}
	case 21:
		return []string{"comunicación", "bienes raíces", "viajes", "vacaciones", "casa", "carga de trabajo", "startups", "finanzas", "sexo", "matrimonio", "concepción", "jardín"}
	case 22:
		return []string{"entrenamiento", "alco", "concepción", "carga", "sueños", "jardín"}
	case 23:
		return []string{"casa", "jardín"}
	case 24:
		return []string{"sexo", "carga", "comunicación", "finanzas", "bienes raíces", "creatividad", "concepción", "hogar", "jardín"}
	case 25:
		return []string{"bienes raíces", "vacaciones", "concepción", "jardín"}
	case 26:
		return []string{"aprendizaje", "matrimonio", "concebir", "hogar", "jardín", "sueños"}
	case 27:
		return []string{"comienzos", "finanzas", "corte de pelo", "finanzas", "viajes", "matrimonio", "jardín"}
	case 28:
		return []string{"creatividad", "concepción", "sueños", "comienzos", "finanzas", "bienes raíces", "hogar", "cargas"}
	case 29:
		return []string{"casa", "corte de pelo"}
	case 30:
		return []string{"descanso", "creatividad", "aprendizaje", "concepción"}
	default:
		panic("moon age error")
	}
}

func (m *Moon) positiveFr() []string {
	switch int(m.Age) {
	case 1:
		return []string{"voyage", "conception"}
	case 2:
		return []string{"débuts", "finance", "apprentissage", "créativité", "voyage", "conception", "maison"}
	case 3:
		return []string{"load", "garden", "sex"}
	case 4:
		return []string{"home", "finance", "vacances", "conception", "jardin"}
	case 5:
		return []string{"alco", "haircut", "load", "garden"}
	case 6:
		return []string{"communication", "apprentissage", "créativité", "maison"}
	case 7:
		return []string{"communication", "rêves", "maison", "repos", "créativité"}
	case 8:
		return []string{"voyage", "immobilier", "vacances", "début", "créativité", "rêves"}
	case 9:
		return []string{"charges", "jardin", "voyage"}
	case 10:
		return []string{"débuts", "vacances", "maison", "finance", "communication", "créativité", "conception", "immobilier", "charge", "jardin"}
	case 11:
		return []string{"voyage", "vacances", "coupe de cheveux", "jardin", "conception", "rêves", "créativité"}
	case 12:
		return []string{"jardin", "rêves", "coupe de cheveux", "socialisation"}
	case 13:
		return []string{"communication", "maison", "sexe", "créativité", "apprentissage", "rêves", "finance", "alcool", "coupe de cheveux"}
	case 14:
		return []string{"débuts", "finance", "créativité", "voyage", "maison", "communication", "apprentissage", "jardin", "charge", "conception"}
	case 15:
		return []string{"maison", "jardin"}
	case 16:
		return []string{"communication", "vacances", "maison", "voyage", "créativité", "mariage", "conception", "rêves"}
	case 17:
		return []string{"communication", "vacances", "sexe", "alco", "créativité", "finance", "immobilier", "voyage", "mariage", "rêves", "coupe de cheveux", "jardin"}
	case 18:
		return []string{"voyages", "créativité", "charges", "rêves"}
	case 19:
		return []string{"voyage", "art", "maison", "jardin", "rêves"}
	case 20:
		return []string{"débuts", "communication", "finance", "immobilier", "créativité", "maison", "voyage", "apprentissage", "jardin"}
	case 21:
		return []string{"communication", "immobilier", "voyage", "vacances", "maison", "charge de travail", "startups", "finance", "sexe", "mariage", "conception", "jardin"}
	case 22:
		return []string{"formation", "alco", "conception", "charge", "rêves", "jardin"}
	case 23:
		return []string{"maison", "jardin"}
	case 24:
		return []string{"sexe", "charge", "communication", "finance", "immobilier", "créativité", "conception", "maison", "jardin"}
	case 25:
		return []string{"immobilier", "vacances", "conception", "jardin"}
	case 26:
		return []string{"apprentissage", "mariage", "conception", "maison", "jardin", "rêves"}
	case 27:
		return []string{"débuts", "finance", "coupe de cheveux", "finance", "voyage", "mariage", "jardin"}
	case 28:
		return []string{"créativité", "conception", "rêves", "débuts", "finance", "immobilier", "maison", "charges"}
	case 29:
		return []string{"maison", "coupe de cheveux"}
	case 30:
		return []string{"repos", "créativité", "apprentissage", "conception"}
	default:
		panic("moon age error")
	}
}

func (m *Moon) positiveIt() []string {
	switch int(m.Age) {
	case 1:
		return []string{"viaggio", "concezione"}
	case 2:
		return []string{"inizio", "finanza", "apprendimento", "creatività", "viaggio", "concepimento", "casa"}
	case 3:
		return []string{"load", "garden", "sex"}
	case 4:
		return []string{"casa", "finanza", "vacanza", "concezione", "giardino"}
	case 5:
		return []string{"alco", "haircut", "load", "garden"}
	case 6:
		return []string{"comunicazione", "apprendimento", "creatività", "casa"}
	case 7:
		return []string{"comunicazione", "sogni", "casa", "riposo", "creatività"}
	case 8:
		return []string{"viaggi", "immobili", "vacanze", "iniziare", "creatività", "sogni"}
	case 9:
		return []string{"loads", "garden", "travel"}
	case 10:
		return []string{"inizio", "vacanza", "casa", "finanza", "comunicazione", "creatività", "concezione", "immobiliare", "carico", "giardino"}
	case 11:
		return []string{"viaggio", "vacanza", "taglio di capelli", "giardino", "concezione", "sogni", "creatività"}
	case 12:
		return []string{"giardino", "sogni", "taglio di capelli", "socializzazione"}
	case 13:
		return []string{"comunicazione", "casa", "sesso", "creatività", "apprendimento", "sogni", "finanza", "alcol", "taglio di capelli"}
	case 14:
		return []string{"inizio", "finanza", "creatività", "viaggio", "casa", "comunicazione", "apprendimento", "giardino", "carico", "concezione"}
	case 15:
		return []string{"casa", "giardino"}
	case 16:
		return []string{"comunicazione", "vacanza", "casa", "viaggio", "creatività", "matrimonio", "concepimento", "sogni"}
	case 17:
		return []string{"comunicazione", "vacanza", "sesso", "alcol", "creatività", "finanza", "immobiliare", "viaggio", "matrimonio", "sogni", "taglio di capelli", "giardino"}
	case 18:
		return []string{"viaggi", "creatività", "carichi", "sogni"}
	case 19:
		return []string{"viaggio", "arte", "casa", "giardino", "sogni"}
	case 20:
		return []string{"inizio", "comunicazione", "finanza", "immobiliare", "creatività", "casa", "viaggio", "apprendimento", "giardino"}
	case 21:
		return []string{"comunicazione", "immobiliare", "viaggio", "vacanza", "casa", "carico di lavoro", "startup", "finanza", "sesso", "matrimonio", "concezione", "giardino"}
	case 22:
		return []string{"training", "alco", "conception", "load", "dreams", "garden"}
	case 23:
		return []string{"casa", "giardino"}
	case 24:
		return []string{"sesso", "carico", "comunicazione", "finanza", "immobiliare", "creatività", "concezione", "casa", "giardino"}
	case 25:
		return []string{"immobile", "vacanza", "concezione", "giardino"}
	case 26:
		return []string{"apprendimento", "matrimonio", "concepimento", "casa", "giardino", "sogni"}
	case 27:
		return []string{"inizio", "finanza", "taglio di capelli", "finanza", "viaggio", "matrimonio", "giardino"}
	case 28:
		return []string{"creatività", "concezione", "sogni", "inizi", "finanza", "immobiliare", "casa", "carichi"}
	case 29:
		return []string{"casa", "taglio di capelli"}
	case 30:
		return []string{"riposo", "creatività", "apprendimento", "concezione"}
	default:
		panic("moon age error")
	}
}

func (m *Moon) positiveEn() []string {
	switch int(m.Age) {
	case 1:
		return []string{"travel", "conception"}
	case 2:
		return []string{"beginnings", "finance", "learning", "creativity", "travel", "conceiving", "home"}
	case 3:
		return []string{"load", "garden", "sex"}
	case 4:
		return []string{"home", "finance", "vacation", "conception", "garden"}
	case 5:
		return []string{"alco", "haircut", "load", "garden"}
	case 6:
		return []string{"communication", "learning", "creativity", "home"}
	case 7:
		return []string{"communication", "dreams", "home", "rest", "creativity"}
	case 8:
		return []string{"travel", "real estate", "vacation", "starting", "creativity", "dreams"}
	case 9:
		return []string{"loads", "garden", "travel"}
	case 10:
		return []string{"beginnings", "vacation", "home", "finance", "communication", "creativity", "conception", "real estate", "load", "garden"}
	case 11:
		return []string{"travel", "vacation", "haircut", "garden", "conception", "dreams", "creativity"}
	case 12:
		return []string{"garden", "dreams", "haircut", "socialization"}
	case 13:
		return []string{"communication", "home", "sex", "creativity", "learning", "dreams", "finance", "alcohol", "haircut"}
	case 14:
		return []string{"beginnings", "finance", "creativity", "travelling", "home", "communication", "learning", "garden", "load", "conception"}
	case 15:
		return []string{"house", "garden"}
	case 16:
		return []string{"communication", "vacation", "home", "travel", "creativity", "marriage", "conception", "dreams"}
	case 17:
		return []string{"communication", "vacation", "sex", "alco", "creativity", "finance", "real estate", "travel", "marriage", "dreams", "haircut", " garden"}
	case 18:
		return []string{"travels", "creativity", "loads", "dreams"}
	case 19:
		return []string{"travel", "art", "home", "garden", "dreams"}
	case 20:
		return []string{"beginnings", "communication", "finance", "real estate", "creativity", "home", "travel", "learning", "garden"}
	case 21:
		return []string{"communication", "real estate", "travel", "vacation", "home", "workload", "startups", "finance", "sex", "marriage", "conception", " garden"}
	case 22:
		return []string{"training", "alco", "conception", "load", "dreams", "garden"}
	case 23:
		return []string{"house", "garden"}
	case 24:
		return []string{"sex", "load", "communication", "finance", "real estate", "creativity", "conception", "home", "garden"}
	case 25:
		return []string{"real estate", "vacation", "conception", "garden"}
	case 26:
		return []string{"learning", "marriage", "conceiving", "home", "garden", "dreams"}
	case 27:
		return []string{"beginnings", "finance", "haircut", "finance", "travel", "marriage", "garden"}
	case 28:
		return []string{"creativity", "conception", "dreams", "beginnings", "finance", "real estate", "home", "loads"}
	case 29:
		return []string{"house", "haircut"}
	case 30:
		return []string{"rest", "creativity", "learning", "conception"}
	default:
		panic("moon age error")
	}
}

func (m *Moon) negativeRu() []string {
	switch int(m.Age) {
	case 1:
		return []string{"начинания", "алко", "споры", "стрижка", "финансы", "общение"}
	case 2:
		return []string{"споры", "сны", "общение", "отдых"}
	case 3:
		return []string{"начинания", "обучение", "финансы", "творчество", "путешествия", "дом", "отдых", "недвижимость"}
	case 4:
		return []string{"общение", "начинания", "алко", "стрижка", "недвижимость"}
	case 5:
		return []string{"начинания", "финансы", "общение", "недвижимость", "путешествия", "обучение", "дом"}
	case 6:
		return []string{"стрижка"}
	case 7:
		return []string{"стрижка"}
	case 8:
		return []string{"алко", "споры", "брак"}
	case 9:
		return []string{"алко", "секс", "сны", "начинания", "сны", "недвижимость"}
	case 10:
		return []string{"путешествия", "споры"}
	case 11:
		return []string{"начинания", "финансы", "общение"}
	case 12:
		return []string{"начинания", "финансы", "алко", "дом", "споры", "путешествия", "недвижимость", "творчество"}
	case 13:
		return []string{"начинания", "дом", "споры", "путешествия", "творчество"}
	case 14:
		return []string{"отдых", "алко", "споры", "сны"}
	case 15:
		return []string{"начинания", "финансы", "алко", "секс", "сны", "споры", "зачатие", "общение", "путешествия", "обучение"}
	case 16:
		return []string{"споры", "секс", "стрижка", "алко"}
	case 17:
		return []string{"общение", "споры", "дом"}
	case 18:
		return []string{"споры", "брак", "начинания", "финансы", "недвижимость", "алко", "секс"}
	case 19:
		return []string{"начинания", "общение", "недвижимость", "споры", "зачатие", "нагрузки"}
	case 20:
		return []string{"стрижка", "отдых", "алко", "зачатие", "сны"}
	case 21:
		return nil
	case 22:
		return []string{"недвижимость", "начинания", "общение", "путешествия"}
	case 23:
		return []string{"общение", "отдых", "секс", "алко", "споры", "брак"}
	case 24:
		return []string{"путешествия", "отдых", "алко", "споры"}
	case 25:
		return []string{"стрижка", "начинания", "общение", "путешествия", "алко", "споры", "брак"}
	case 26:
		return []string{"начинания", "общение", "недвижимость", "алко", "споры", "путешествия", "отдых", "секс", "стрижка"}
	case 27:
		return []string{"алко", "сны"}
	case 28:
		return []string{"стрижка", "алко"}
	case 29:
		return []string{"начинания", "общение", "финансы", "недвижимость", "алко", "споры", "секс", "зачатие", "сны"}
	case 30:
		return []string{"стрижка", "общение", "путешествия", "алко", "споры", "сад"}
	default:
		panic("moon age error")
	}
}

func (m *Moon) negativeDe() []string {
	switch int(m.Age) {
	case 1:
		return []string{"Anfänge", "Alco", "Streitigkeiten", "Haarschnitt", "Finanzen", "Kommunikation"}
	case 2:
		return []string{"argumente", "träume", "kommunikation", "ruhe"}
	case 3:
		return []string{"Anfänge", "Lernen", "Finanzen", "Kreativität", "Reisen", "Zuhause", "Urlaub", "Immobilien"}
	case 4:
		return []string{"Kommunikation", "Anfänge", "Alco", "Haarschnitt", "Immobilien"}
	case 5:
		return []string{"Start", "Finanzen", "Kommunikation", "Immobilien", "Reisen", "Training", "Zuhause"}
	case 6:
		return []string{"haircut"}
	case 7:
		return []string{"haircut"}
	case 8:
		return []string{"alkohol", "argumente", "ehe"}
	case 9:
		return []string{"Alco", "Sex", "Träume", "Anfänge", "Träume", "Immobilien"}
	case 10:
		return []string{"Reisen", "Streitigkeiten"}
	case 11:
		return []string{"Anfänge", "Finanzen", "Kommunikation"}
	case 12:
		return []string{"Anfänge", "Finanzen", "Alco", "Zuhause", "Streitigkeiten", "Reisen", "Immobilien", "Kreativität"}
	case 13:
		return []string{"Anfänge", "Heimat", "Streitigkeiten", "Reisen", "Kreativität"}
	case 14:
		return []string{"rest", "alco", "arguments", "dreams"}
	case 15:
		return []string{"Anfänge", "Finanzen", "Alco", "Sex", "Träume", "Argumente", "Konzeption", "Kommunikation", "Reisen", "Lernen"}
	case 16:
		return []string{"argument", "sex", "haircut", "alco"}
	case 17:
		return []string{"Kommunikation", "Argumente", "Haus"}
	case 18:
		return []string{"Streitigkeiten", "Ehe", "Unternehmungen", "Finanzen", "Immobilien", "Alkohol", "Sex"}
	case 19:
		return []string{"Anfänge", "Kommunikation", "Immobilien", "Streitigkeiten", "Konzeption", "Lasten"}
	case 20:
		return []string{"haircut", "rest", "alco", "konzeption", "träume"}
	case 21:
		return nil
	case 22:
		return []string{"Immobilien", "Anfänge", "Kommunikation", "Reisen"}
	case 23:
		return []string{"Kommunikation", "Urlaub", "Sex", "Alkohol", "Argumente", "Ehe"}
	case 24:
		return []string{"travel", "vacation", "alco", "spores"}
	case 25:
		return []string{"Haarschnitt", "Anfänge", "Kommunikation", "Reisen", "Alkohol", "Argumente", "Ehe"}
	case 26:
		return []string{"Anfänge", "Kommunikation", "Immobilien", "Alkohol", "Argumente", "Reisen", "Urlaub", "Sex", "Haarschnitt"}
	case 27:
		return []string{"alco", "dreams"}
	case 28:
		return []string{"haircut", "alco"}
	case 29:
		return []string{"Anfänge", "Kommunikation", "Finanzen", "Immobilien", "Alkohol", "Argumente", "Sex", "Empfängnis", "Träume"}
	case 30:
		return []string{"haircut", "social", "travel", "alco", "arguments", "garden"}
	default:
		panic("moon age error")
	}
}

func (m *Moon) negativeEs() []string {
	switch int(m.Age) {
	case 1:
		return []string{"comienzos", "alco", "disputas", "corte de pelo", "finanzas", "comunicación"}
	case 2:
		return []string{"argumentos", "sueños", "comunicación", "descanso"}
	case 3:
		return []string{"comienzos", "aprendizaje", "finanzas", "creatividad", "viajes", "casa", "vacaciones", "bienes raíces"}
	case 4:
		return []string{"comunicación", "comienzos", "alco", "corte de pelo", "bienes raíces"}
	case 5:
		return []string{"inicio", "finanzas", "comunicación", "bienes raíces", "viajes", "capacitación", "casa"}
	case 6:
		return []string{"corte de pelo"}
	case 7:
		return []string{"corte de pelo"}
	case 8:
		return []string{"alcohol", "argumentos", "matrimonio"}
	case 9:
		return []string{"alco", "sexo", "sueños", "comienzos", "sueños", "bienes raíces"}
	case 10:
		return []string{"viajes", "disputas"}
	case 11:
		return []string{"comienzos", "finanzas", "comunicación"}
	case 12:
		return []string{"comienzos", "finanzas", "alco", "casa", "disputas", "viajes", "bienes raíces", "creatividad"}
	case 13:
		return []string{"comienzos", "casa", "disputas", "viajes", "creatividad"}
	case 14:
		return []string{"descanso", "alco", "argumentos", "sueños"}
	case 15:
		return []string{"comienzos", "finanzas", "alco", "sexo", "sueños", "argumentos", "concepción", "comunicación", "viaje", "aprendizaje"}
	case 16:
		return []string{"argumento", "sexo", "corte de pelo", "alcohólico"}
	case 17:
		return []string{"comunicación", "argumentos", "casa"}
	case 18:
		return []string{"disputas", "matrimonio", "compromisos", "finanzas", "bienes raíces", "alcohol", "sexo"}
	case 19:
		return []string{"comienzos", "comunicación", "bienes raíces", "disputas", "concepción", "cargas"}
	case 20:
		return []string{"corte de pelo", "descanso", "alcohólico", "concepción", "sueños"}
	case 21:
		return nil
	case 22:
		return []string{"bienes raíces", "comienzos", "comunicación", "viajes"}
	case 23:
		return []string{"comunicación", "vacaciones", "sexo", "alcohol", "argumentos", "matrimonio"}
	case 24:
		return []string{"viaje", "vacaciones", "alcohólico", "esporas"}
	case 25:
		return []string{"corte de pelo", "comienzos", "comunicación", "viajes", "alcohol", "argumentos", "matrimonio"}
	case 26:
		return []string{"comienzos", "comunicación", "bienes raíces", "alcohol", "argumentos", "viajes", "vacaciones", "sexo", "corte de pelo"}
	case 27:
		return []string{"alco", "sueños"}
	case 28:
		return []string{"corte de pelo", "alcohólico"}
	case 29:
		return []string{"comienzos", "comunicación", "finanzas", "bienes raíces", "alcohol", "argumentos", "sexo", "concepción", "sueños"}
	case 30:
		return []string{"corte de pelo", "social", "viaje", "alcohólico", "argumentos", "jardín"}
	default:
		panic("moon age error")
	}
}

func (m *Moon) negativeFr() []string {
	switch int(m.Age) {
	case 1:
		return []string{"débuts", "alco", "litiges", "coupe de cheveux", "finance", "communication"}
	case 2:
		return []string{"arguments", "rêves", "communication", "repos"}
	case 3:
		return []string{"débuts", "apprentissage", "finance", "créativité", "voyage", "maison", "vacances", "immobilier"}
	case 4:
		return []string{"communication", "débuts", "alco", "coupe de cheveux", "immobilier"}
	case 5:
		return []string{"starting", "finance", "communication", "immobilier", "voyage", "formation", "maison"}
	case 6:
		return []string{"coupe de cheveux"}
	case 7:
		return []string{"coupe de cheveux"}
	case 8:
		return []string{"alcool", "arguments", "mariage"}
	case 9:
		return []string{"alco", "sexe", "rêves", "débuts", "rêves", "immobilier"}
	case 10:
		return []string{"voyages", "contestations"}
	case 11:
		return []string{"débuts", "finance", "communication"}
	case 12:
		return []string{"débuts", "finance", "alco", "maison", "litiges", "voyages", "immobilier", "créativité"}
	case 13:
		return []string{"débuts", "maison", "contestations", "voyages", "créativité"}
	case 14:
		return []string{"rest", "alco", "arguments", "dreams"}
	case 15:
		return []string{"débuts", "finance", "alco", "sexe", "rêves", "arguments", "conception", "communication", "voyage", "apprentissage"}
	case 16:
		return []string{"argument", "sexe", "coupe de cheveux", "alco"}
	case 17:
		return []string{"communication", "arguments", "maison"}
	case 18:
		return []string{"litiges", "mariage", "engagements", "finance", "immobilier", "alcool", "sexe"}
	case 19:
		return []string{"débuts", "communication", "immobilier", "litiges", "conception", "charges"}
	case 20:
		return []string{"coupe de cheveux", "repos", "alco", "conception", "rêves"}
	case 21:
		return nil
	case 22:
		return []string{"immobilier", "débuts", "communication", "voyage"}
	case 23:
		return []string{"communication", "vacances", "sexe", "alcool", "arguments", "mariage"}
	case 24:
		return []string{"voyage", "vacances", "alco", "spores"}
	case 25:
		return []string{"coupe de cheveux", "débuts", "communication", "voyage", "alcool", "disputes", "mariage"}
	case 26:
		return []string{"débuts", "communication", "immobilier", "alcool", "disputes", "voyage", "vacances", "sexe", "coupe de cheveux"}
	case 27:
		return []string{"alco", "rêves"}
	case 28:
		return []string{"coupe de cheveux", "alco"}
	case 29:
		return []string{"débuts", "communication", "finance", "immobilier", "alcool", "disputes", "sexe", "conception", "rêves"}
	case 30:
		return []string{"haircut", "social", "travel", "alco", "arguments", "garden"}
	default:
		panic("moon age error")
	}
}

func (m *Moon) negativeIt() []string {
	switch int(m.Age) {
	case 1:
		return []string{"inizi", "alco", "controversie", "taglio di capelli", "finanza", "comunicazione"}
	case 2:
		return []string{"argomenti", "sogni", "comunicazione", "riposo"}
	case 3:
		return []string{"inizio", "apprendimento", "finanza", "creatività", "viaggio", "casa", "vacanza", "immobiliare"}
	case 4:
		return []string{"communication", "beginnings", "alco", "haircut", "real estate"}
	case 5:
		return []string{"starting", "finance", "communication", "real estate", "travel", "training", "home"}
	case 6:
		return []string{"taglio di capelli"}
	case 7:
		return []string{"taglio di capelli"}
	case 8:
		return []string{"alcol", "argomenti", "matrimonio"}
	case 9:
		return []string{"alco", "sex", "dreams", "beginnings", "dreams", "real estate"}
	case 10:
		return []string{"viaggi", "controversie"}
	case 11:
		return []string{"inizio", "finanza", "comunicazione"}
	case 12:
		return []string{"inizi", "finanza", "alco", "casa", "controversie", "viaggi", "immobiliare", "creatività"}
	case 13:
		return []string{"inizi", "home", "controversie", "viaggi", "creatività"}
	case 14:
		return []string{"riposo", "alco", "argomenti", "sogni"}
	case 15:
		return []string{"inizi", "finanza", "alco", "sesso", "sogni", "argomenti", "concezione", "comunicazione", "viaggio", "apprendimento"}
	case 16:
		return []string{"argomento", "sesso", "taglio di capelli", "alco"}
	case 17:
		return []string{"comunicazione", "argomenti", "casa"}
	case 18:
		return []string{"controversie", "matrimonio", "impegni", "finanza", "immobili", "alcol", "sesso"}
	case 19:
		return []string{"inizi", "comunicazione", "real estate", "controversie", "concezione", "carichi"}
	case 20:
		return []string{"haircut", "rest", "alco", "conception", "dreams"}
	case 21:
		return nil
	case 22:
		return []string{"immobiliare", "inizio", "comunicazione", "viaggio"}
	case 23:
		return []string{"comunicazione", "vacanza", "sesso", "alcol", "argomenti", "matrimonio"}
	case 24:
		return []string{"viaggio", "vacanza", "alco", "spore"}
	case 25:
		return []string{"haircut", "beginnings", "communication", "travel", "alcohol", "litigi", "matrimonio"}
	case 26:
		return []string{"inizio", "comunicazione", "immobiliare", "alcol", "argomenti", "viaggio", "vacanza", "sesso", "taglio di capelli"}
	case 27:
		return []string{"alco", "sogni"}
	case 28:
		return []string{"haircut", "alco"}
	case 29:
		return []string{"inizio", "comunicazione", "finanza", "immobiliare", "alcol", "argomenti", "sesso", "concezione", "sogni"}
	case 30:
		return []string{"haircut", "social", "travel", "alco", "arguments", "garden"}
	default:
		panic("moon age error")
	}
}

func (m *Moon) negativeUa() []string {
	switch int(m.Age) {
	case 1:
		return []string{"починання", "алко", "спори", "стрижка", "фінанси", "спілкування"}
	case 2:
		return []string{"суперечки", "сни", "спілкування", "відпочинок"}
	case 3:
		return []string{"починання", "навчання", "фінанси", "творчість", "подорожі", "будинок", "відпочинок", "нерухомість"}
	case 4:
		return []string{"спілкування", "починання", "алко", "стрижка", "нерухомість"}
	case 5:
		return []string{"починання", "фінанси", "спілкування", "нерухомість", "подорожі", "навчання", "будинок"}
	case 6:
		return []string{"стрижка"}
	case 7:
		return []string{"стрижка"}
	case 8:
		return []string{"алко", "спори", "шлюб"}
	case 9:
		return []string{"алко", "секс", "сни", "починання", "сни", "нерухомість"}
	case 10:
		return []string{"подорожі", "суперечки"}
	case 11:
		return []string{"починання", "фінанси", "спілкування"}
	case 12:
		return []string{"починання", "фінанси", "алко", "будинок", "спори", "подорожі", "нерухомість", "творчість"}
	case 13:
		return []string{"починання", "будинок", "суперечки", "подорожі", "творчість"}
	case 14:
		return []string{"відпочинок", "алко", "спори", "сни"}
	case 15:
		return []string{"починання", "фінанси", "алко", "секс", "сни", "суперечки", "зачаття", "спілкування", "подорожі", "навчання"}
	case 16:
		return []string{"спори", "секс", "стрижка", "алко"}
	case 17:
		return []string{"спілкування", "спори", "дім"}
	case 18:
		return []string{"спори", "шлюб", "починання", "фінанси", "нерухомість", "алко", "секс"}
	case 19:
		return []string{"починання", "спілкування", "нерухомість", "спори", "зачаття", "навантаження"}
	case 20:
		return []string{"стрижка", "відпочинок", "алко", "зачаття", "сни"}
	case 21:
		return nil
	case 22:
		return []string{"нерухомість", "починання", "спілкування", "подорожі"}
	case 23:
		return []string{"спілкування", "відпочинок", "секс", "алко", "спори", "шлюб"}
	case 24:
		return []string{"подорожі", "відпочинок", "алко", "суперечки"}
	case 25:
		return []string{"стрижка", "починання", "спілкування", "подорожі", "алко", "спори", "шлюб"}
	case 26:
		return []string{"починання", "спілкування", "нерухомість", "алко", "спори", "подорожі", "відпочинок", "секс", "стрижка"}
	case 27:
		return []string{"алко", "сни"}
	case 28:
		return []string{"стрижка", "алко"}
	case 29:
		return []string{"починання", "спілкування", "фінанси", "нерухомість", "алко", "спори", "секс", "зачаття", "сни"}
	case 30:
		return []string{"стрижка", "спілкування", "подорожі", "алко", "спори", "сад"}
	default:
		panic("moon age error")
	}
}

func (m *Moon) negativeEn() []string {
	switch int(m.Age) {
	case 1:
		return []string{"beginnings", "alco", "disputes", "haircut", "finance", "communication"}
	case 2:
		return []string{"arguments", "dreams", "communication", "rest"}
	case 3:
		return []string{"beginnings", "learning", "finance", "creativity", "travel", "home", "vacation", "real estate"}
	case 4:
		return []string{"communication", "beginnings", "alco", "haircut", "real estate"}
	case 5:
		return []string{"starting", "finance", "communication", "real estate", "travel", "training", "home"}
	case 6:
		return []string{"haircut"}
	case 7:
		return []string{"haircut"}
	case 8:
		return []string{"alcohol", "arguments", "marriage"}
	case 9:
		return []string{"alco", "sex", "dreams", "beginnings", "dreams", "real estate"}
	case 10:
		return []string{"travels", "disputes"}
	case 11:
		return []string{"beginnings", "finance", "communication"}
	case 12:
		return []string{"beginnings", "finance", "alco", "home", "disputes", "travel", "real estate", "creativity"}
	case 13:
		return []string{"beginnings", "home", "disputes", "travel", "creativity"}
	case 14:
		return []string{"rest", "alco", "arguments", "dreams"}
	case 15:
		return []string{"beginnings", "finance", "alco", "sex", "dreams", "arguments", "conception", "communication", "travel", "learning"}
	case 16:
		return []string{"argument", "sex", "haircut", "alco"}
	case 17:
		return []string{"communication", "arguments", "house"}
	case 18:
		return []string{"disputes", "marriage", "undertakings", "finance", "real estate", "alcohol", "sex"}
	case 19:
		return []string{"beginnings", "communication", "real estate", "disputes", "conception", "loads"}
	case 20:
		return []string{"haircut", "rest", "alco", "conception", "dreams"}
	case 21:
		return nil
	case 22:
		return []string{"real estate", "beginnings", "communication", "travel"}
	case 23:
		return []string{"communication", "vacation", "sex", "alcohol", "arguments", "marriage"}
	case 24:
		return []string{"travel", "vacation", "alco", "spores"}
	case 25:
		return []string{"haircut", "beginnings", "communication", "travel", "alcohol", "arguments", "marriage"}
	case 26:
		return []string{"beginnings", "communication", "real estate", "alcohol", "arguments", "travel", "vacation", "sex", "haircut"}
	case 27:
		return []string{"alco", "dreams"}
	case 28:
		return []string{"haircut", "alco"}
	case 29:
		return []string{"beginnings", "communication", "finance", "real estate", "alcohol", "arguments", "sex", "conception", "dreams"}
	case 30:
		return []string{"haircut", "social", "travel", "alco", "arguments", "garden"}
	default:
		panic("moon age error")
	}
}
