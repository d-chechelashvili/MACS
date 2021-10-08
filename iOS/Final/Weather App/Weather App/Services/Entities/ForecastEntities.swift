//
//  ForecastEntities.swift
//  Weather App
//
//  Created by dato che on 2/7/21.
//

import Foundation


struct FiveDayForecast: Codable{
    let list: [TimeStamp]
}

struct TimeStamp: Codable{
    let main: MainInfo
    let weather: [WeatherInfo]
    let dt_txt: String
}

