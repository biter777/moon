// Package moon - package moon
package moon

import (
	"math"
)

const (
	// Density - The average density of the moon, grams per cubic centimeter.
	Density = 3.3464
	// AverageDiameter - Average diameter, meters
	AverageDiameter = 3474200.0
	// AverageDistance - Average distance from earth to moon, meters
	AverageDistance = 384401000.0
	// Gravity - Acceleration of gravity, meters per second squared
	Gravity = 1.62
	// Mass - The mass of the moon, kilograms
	Mass = 7.3477e+22
	// AverageRadius - Average radius, meters
	AverageRadius = AverageDiameter / 2
	// Square - Square, meters
	Square = 4 * math.Pi * AverageRadius * AverageRadius
	// SynodicMonth - Synodic month (new Moon to new Moon), days (1 day = 24h)
	SynodicMonth = 29.53058868
)

const (
	// Astronomical constants
	epoch = 2444238.5 // 1989 January 0.0
	// Constants defining the Sun's apparent orbit
	elonge    = 278.833540 // Ecliptic longitude of the Sun at epoch 1980.0
	elongp    = 282.596403 // Ecliptic longitude of the Sun at perigee
	eccent    = 0.016718   // Eccentricity of Earth's orbit
	sunsmax   = 1.495985e8 // Sun's angular size, degrees, at semi-major axis distance
	sunangsiz = 0.533128
	// Elements of the Moon's orbit, epoch 1980.0
	mmlong  = 64.975464              // Moon's mean longitude at the epoch
	mmlongp = 349.383063             // Mean longitude of the perigee at the epoch
	mecc    = 0.054900               // Eccentricity of the Moon's orbit
	mangsiz = 0.5181                 // Moon's angular size at distance a from Earth
	msmax   = AverageDistance / 1000 // Semi-major axis of Moon's orbit in km
)
