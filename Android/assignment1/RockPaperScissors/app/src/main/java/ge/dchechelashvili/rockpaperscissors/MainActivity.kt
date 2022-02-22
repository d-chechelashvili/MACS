package ge.dchechelashvili.rockpaperscissors

import android.os.Bundle
import android.view.View
import android.widget.ImageView
import android.widget.TextView
import androidx.appcompat.app.AppCompatActivity
import androidx.core.content.ContextCompat
import kotlin.random.Random

enum class Choice {
    ROCK,
    PAPER,
    SCISSORS
}

enum class Result {
    COMPUTER,
    PLAYER,
    DRAW
}

class MainActivity : AppCompatActivity() {
    private var gameStarted = false
    private var playerScore = 0
    private var computerScore = 0

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activity_main)
    }

    fun buttonClicked(view: View) {
        val computerChoice = Choice.values()[Random.nextInt(3)]
        val playerChoice = view.id
        val (result, computerImage, playerImage) = getOutcome(computerChoice, playerChoice)

        updateScore(result)
        updateImages(computerImage, playerImage)

        if (!gameStarted) {
            gameStarted = true
            findViewById<TextView>(R.id.startTextView).visibility = View.INVISIBLE
            findViewById<ImageView>(R.id.computerImage).visibility = View.VISIBLE
            findViewById<ImageView>(R.id.playerImage).visibility = View.VISIBLE
        }
    }

    private fun updateImages(computerImage: Int, playerImage: Int) {
        findViewById<ImageView>(R.id.playerImage).setImageResource(playerImage)
        findViewById<ImageView>(R.id.computerImage).setImageResource(computerImage)
    }


    private fun updateScore(result: Result) {
        val drawColor = ContextCompat.getColor(applicationContext, R.color.drawColor)
        val winColor = ContextCompat.getColor(applicationContext, R.color.winColor)
        val loseColor = ContextCompat.getColor(applicationContext, R.color.loseColor)

        var computerColor = drawColor
        var playerColor = drawColor

        when (result) {
            Result.COMPUTER -> {
                computerColor = winColor
                playerColor = loseColor
                computerScore++
            }
            Result.PLAYER -> {
                computerColor = loseColor
                playerColor = winColor
                playerScore++
            }
            else -> {

            }
        }

        findViewById<TextView>(R.id.playerScore).setTextColor(playerColor)
        findViewById<TextView>(R.id.playerScore).text = playerScore.toString()
        findViewById<TextView>(R.id.computerScore).setTextColor(computerColor)
        findViewById<TextView>(R.id.computerScore).text = computerScore.toString()
    }

    private fun getOutcome(computerChoice: Choice, playerChoice: Int): Triple<Result, Int, Int> {
        var computerImage = R.drawable.rock
        var playerImage = R.drawable.rock
        var result = Result.COMPUTER

        if (playerChoice == R.id.rockButton) {
            if (computerChoice == Choice.ROCK) {
                result = Result.DRAW
            } else if (computerChoice == Choice.PAPER) {
                computerImage = R.drawable.paper
            } else {
                computerImage = R.drawable.scissors
                result = Result.PLAYER
            }
        } else if (playerChoice == R.id.paperButton) {
            playerImage = R.drawable.paper
            if (computerChoice == Choice.ROCK) {
                result = Result.PLAYER
            } else if (computerChoice == Choice.PAPER) {
                computerImage = R.drawable.paper
                result = Result.DRAW
            } else {
                computerImage = R.drawable.scissors
            }
        } else if (playerChoice == R.id.scisssorButton) {
            playerImage = R.drawable.scissors
            if (computerChoice == Choice.ROCK) {
            } else if (computerChoice == Choice.PAPER) {
                computerImage = R.drawable.paper
                result = Result.PLAYER
            } else {
                computerImage = R.drawable.scissors
                result = Result.DRAW
            }
        }
        return Triple(result, computerImage, playerImage)
    }
}