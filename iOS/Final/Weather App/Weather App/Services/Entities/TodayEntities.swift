//
//  TodayEntities.swift
//  Weather App
//
//  Created by dato che on 2/5/21.
//

import Foundation

struct ServiceResponse: Codable{
    let weather: [WeatherInfo]
    let main: MainInfo
    let wind: WindInfo
    let clouds: CloudInfo
    let sys: SysInfo
    let name: String
}

struct WeatherInfo: Codable{
    let main: String
    let description: String
    let icon: String
}

struct MainInfo: Codable{
    let temp: Double
    let humidity: Int
}

struct WindInfo: Codable{
    let speed: Double
    let deg: Int
}

struct CloudInfo: Codable{
    let all: Int
}

struct SysInfo: Codable{
    let country: String
}
